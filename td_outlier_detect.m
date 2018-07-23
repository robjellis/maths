function [Xmat stats] = td_outlier_detect(L, R, ex_cap, ex_method, ex_inds_import, del_itis, long_run_method, loc_meth, outlier_meth, outlier_thr, string_type, run_dur_thr, gap_dur_thr, do_slope, plot_it,fnum,fname)

% td_outlier_detect(L,R,tau,ex_cap,run_dur_thr,gap_dur_thr)
%
% takes in vectors of actual *timestamps* (continuous from 0, in seconds)
% or *displacements* (continuous from starting location, in cm)
%
% if R is empty, then we assume we just do calcs on the first vector of
% timestamps (e.g., can work for HR data, etc)

% "loc_meth" will either be 'kde_botev', 'kde_matlab' (for BEATS), 'tm' for trimean (for Physionet)
%
% "outlier_meth" is 'pc' for percent change, 'ad' for absolute difference, etc.
%
% "ex_method"
% 1 = Hausdorff auto
% 2 = Hausdorff + manual
% 3 = RJE auto
% 4 = RJE + manual
% 5 = no exclusions

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% perform the location comparison?
do_loc_compare = 0;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(R)
    num_chan = 1;
    L = L(:);
    all_ts = L;
    
    sort_ts = sort(all_ts); % just to be safe
    sort_chan = ones(size(L));
else
    
    num_chan = 2;
    R = R(:);
    all_ts = [L; R]; % unsorted (leave like this)
    all_chan = [ones(size(L)); ones(size(R))+1];
    Rdiff = diff(R);
    Rdiff = [NaN; Rdiff]; % used for plotting
    
    % get sorted versions too
    tt = [all_ts all_chan];
    tt = sortrows(tt,1);
    sort_ts = tt(:,1);
    sort_chan = tt(:,2);
end
    
% keep this here
Ldiff = diff(L);
Ldiff = [NaN; Ldiff]; % used for plotting
Mdiff = diff(sort_ts);
   
if nargin < 3
    ex_cap = 0;
    ex_method = 3;
    ex_inds_import = 0;
    del_itis = 0;
    long_run_method = 'c'; % concatenate after removing outliers
    loc_meth = 'tm'; % trimean
    outlier_meth = 'pc';
    outlier_thr = 50;
    string_type = 's'; % string of events defined either in 's' seconds or 'e' events
    run_dur_thr = 10; % at least N continuous seconds/events 
    gap_dur_thr = Inf; % any duration gap is fine
    do_slope = 0;
    plot_it = 0;
    fnum = 300;
    fname = 'File';
end

%% Hausdorff method: +/- 3 SD away from the median

% stdev
Mdiff_med = median(Mdiff); % note: all values are still included (large and small outliers)
Mdiff_std = std(Mdiff);

% traditional method
loH = Mdiff > Mdiff_med - 3*Mdiff_std;
hiH = Mdiff < Mdiff_med + 3*Mdiff_std;

H_val = min(loH,hiH);
H_val = [NaN; H_val]; % add the missing case *after* since we don't want it to affect median or mad
% note: this doesn't take alternating steps into account


%% rje method

% 0. merge the series

    % performed above

% 2. find the location of the entire series

if do_loc_compare == 1
    loc_compare(Mdiff,1);
end

if strcmp(loc_meth,'bf')
    % brute force method: find the value from which the most number of values have a percentage deviation of <= thr

    brute_plot = 1; % 1 means plot 
    brute_thr = 1; % need this small
    min_prc = 0; % simplest
    brute_out = cent_tend_brute(Mdiff,brute_thr,min_prc,brute_plot);
    
    loc = brute_out.brute_location;

elseif strcmp(loc_meth,'kde_matlab')

    % simple rule: if we don't have enough unique, then just use mode
    prc_unique = 100 * numel(unique(Mdiff)) / numel(Mdiff);
    
    if prc_unique < 50
        % just use mode to be safe
        loc = mode(Mdiff);
    else
        npoints = 2^10;

        range = max(Mdiff) - min(Mdiff);
        MIN = min(Mdiff) - range/4; 
        MAX = max(Mdiff) + range/4;

        [f,xi] = ksdensity(Mdiff,linspace(MIN,MAX,npoints)); 

        loc = min(xi(f == max(f))); % take min so we get the faster rate
    end
    
elseif strcmp(loc_meth,'kde_botev')

    npoints = 2^10; % default is 2^12 
    
    [bandwidth,density,xmesh] = kde(Mdiff,npoints);
    
    loc = min(xmesh(density == max(density))); % take min so we get the lowest peak

elseif strcmp(loc_meth,'tm')
    % http://en.wikipedia.org/wiki/Trimean
    
    % 2b. for Physionet, we use the trimean rather than the median, since
    % it will better reflect the positive skew of the distribution
    
    % we can even do a manual scrub of obvious outliers to give an even
    % better estimate of the location. Use the "4x ratio relative to median rule" here
    
%     med_Mdiff = median(Mdiff);
%     
%     Mdiff_scrub = Mdiff(Mdiff             > med_Mdiff * 4);
%     Mdiff_scrub = Mdiff_scrub(Mdiff_scrub < med_Mdiff / 4);
    
    loc = (prctile_nist(Mdiff,25) + 2 * prctile_nist(Mdiff,50) + prctile_nist(Mdiff,75)) / 4; % Tukey formula

elseif strcmp(loc_meth,'mean')
    loc = mean(Mdiff);
end

% important: adjust the indices of Mdiff to match original event series;
Mdiff = [loc; Mdiff];

%% now we choose the desired outlier thresholding method
% 2013.10.11 - modified now such that it is a TWO-STEP PROCESS

% % first we get the "grid" - the theoretical timestamps based on loc
% grid = repmat(loc,num_ts,1);
% grid = cumsum(grid); % timestamps
% grid_ts = [0; grid]; % first timestamp is zero
% 
% % now we need to get the grid aligned with the actual timestamps in terms of signed error
% err = sort_ts - grid_ts;
% grid_ts = grid_ts + err;

do_suc = 1; % default, leave this here

% adjust the indices of Mdiff

if strcmp(outlier_meth,'ad')  
    % need to adjust the indices

    loc_trans = abs(Mdiff);
    do_suc = 0;
    
    stable_bin = loc_trans < outlier_thr;
elseif strcmp(outlier_meth,'gt')
    % all first order diffs must be greater than some value; useful for
    % rapid tapping outlier detection
    loc_trans = [];
    do_suc = 0;
    
    stable_bin = Mdiff > outlier_thr; % to find successive timestamps that are close in time
    
elseif strcmp(outlier_meth,'mult_loc')
    
    % a multiple of the location
    loc_trans = abs(Mdiff);
    do_suc = 0;
    
    thr = outlier_thr * loc;
    stable_bin = loc_trans < thr;   

elseif strcmp(outlier_meth,'add_loc')
    
    % add value to location
    
elseif strcmp(outlier_meth,'ld')
    do_suc = 0;
    loc_trans = abs(Mdiff);
    do_suc = 0;
    
elseif strcmp(outlier_meth,'pc') || strcmp(outlier_meth,'pc_suc')
    do_suc = 1;
    % revised 2014.07.19 so that we can opt to just have the successive transformation
    if strcmp(outlier_meth,'pc')
    % 1. deviation from location
        loc_trans = 100 * abs(Mdiff - loc) / loc;

        loc_trans = loc_trans(:);
        %loc_trans = [NaN; loc_trans]; % so index locations are correct (for timestamps)
        stable_bin = loc_trans <= outlier_thr; % stable_bin(1) will always be 0 (excluded)
    else
        stable_bin = ones(size(Mdiff));
        %stable_bin = [1; stable_bin]; % so index locations are correct
        loc_trans = nan(size(stable_bin)); % placeholder
    end
    
    % 2. successive 
    suc_trans_all =  nan(numel(Mdiff-1),1);
    suc_trans_good = nan(numel(Mdiff-1),1);
    
    for i = 1:numel(Mdiff)-1
        if stable_bin(i) == 1 && stable_bin(i+1) == 1
            %suc_trans = 100 * abs(Mdiff(i+1) - Mdiff(i)) / (Mdiff(i+1) + Mdiff(i));
            suc_trans = 100 * abs(Mdiff(i+1) - Mdiff(i)) / Mdiff(i);
            suc_trans_all(i) = suc_trans;
            if suc_trans <= outlier_thr
               % OK; stable_bin stays at 1
               
               % store this value
               suc_trans_good(i) = suc_trans; % we only store "good" values
            else
               stable_bin(i+1) = 0;
            end
        else
            % do nothing
        end
    

    end 
    % now simply remove NaN from suc_trans_good
    suc_trans_good(isnan(suc_trans_good)) = [];
    
elseif strcmp(outlier_meth,'pd') % percent difference
    
    % 1. deviation from location
    loc_trans = 100 * abs(Mdiff - loc) / (0.5 * (Mdiff + loc));

    loc_trans = loc_trans(:);
    loc_trans = [NaN; loc_trans]; % so index locations are correct (for timestamps)
    stable_bin = loc_trans <= outlier_thr; % stable_bin(1) will always be 0 (excluded)

    % 2. successive 
    suc_trans_all = nan(numel(Mdiff-1),1);
    suc_trans_good = nan(numel(Mdiff-1),1);

    for i = 1:numel(Mdiff)-1
        if stable_bin(i) == 1 && stable_bin(i+1) == 1
            suc_trans = 100 * abs(Mdiff(i+1) - Mdiff(i)) / (0.5 * (Mdiff(i+1) + Mdiff(i)));
            suc_trans_all(i) = suc_trans;
            if suc_trans <= outlier_thr
               % OK; stable_bin stays at 1
               suc_trans_good(i) = suc_trans; % we only store "good" values
            else
               stable_bin(i+1) = 0;
            end
        else
            % do nothing
        end
    end        
    
    % now simply remove NaN from suc_trans_good
    suc_trans_good(isnan(suc_trans_good)) = [];
    
end 

    

% 6. the final included (= 1) and excluded (= 0) cases

E_val = stable_bin;

% now find loc_transf_max and suc_transf_max 

if do_suc == 1
    loc_transf_max = max(loc_trans(stable_bin == 1));
    suc_transf_max = max(suc_trans_good);
else
    loc_transf_max = [];
    suc_transf_max = [];
end

%% pool the methods

% note: H_val and E_val must be 0 or 1 for all values for long_run_ts.m to work
Xmat(:,[1 2 3 4]) = [sort_ts sort_chan H_val E_val]; % four columns

% cut out the outer cases for simplicity

beg_ind = ex_cap+1;
end_ind = numel(all_ts)-ex_cap;

Xmat = Xmat(beg_ind:end_ind,:); % cut these cases out

% leave this here
mT = Xmat(:,1);
dmT = diff(mT);  % inter-event intervals
dmT = [NaN; dmT]; % to match length of mT



%% show plots
if plot_it >= 1 % keep this like this
    
    if ex_method == 2 || ex_method == 4
        real_time = 0; % need actual indicies, not real time
    else
        real_time = 1; % plot in seconds (default)
    end
 
    % force it
    real_time = 0;
    
    figure(50)
    nplots = 4;
    
    % plot L and R directly
    
    subplot(nplots,1,1)
%     if real_time == 0
%         plot(Ldiff, 'b','Marker','.')
%         hold on 
%         if num_chan == 2
%             plot(Rdiff, 'r','Marker','.')
%         end
%         xlabel('Ordinal strides')
%     elseif real_time == 1
        plot(L,Ldiff, 'b','Marker','.')
        hold on
        if num_chan == 2
            plot(R,Rdiff, 'r','Marker','.')
        end
        xlabel('Cumulative time (s) or displacement (cm)')
        title('[this figure from td_outlier_detect.m]','Interpreter','none')
%     end
    
    ylabel('Inter-STRIDE Int.')
    hold off
    
    if size(fname) > 0
        title(['File ' num2str(fnum) ': ' num2str(fname)],'Interpreter','none')
    end
    
    
    subplot(nplots,1,2)

    if real_time == 0
        plot(dmT,'k','Marker','.')
        hold on
        plot([1 numel(dmT)],[loc loc],'c') % location
        
    elseif real_time == 1
        plot(mT,dmT,'k','Marker','.')
        hold on
        plot([min(dmT) max(dmT)],[loc loc],'c') % location
    end
    
    % get the H outliers: col3
    H = Xmat(:,3);
    xH =  mT;
    yH = dmT;
    
    xH(H == 1) = NaN; 
    yH(H == 1) = NaN; 

    % get the E outliers: col4
    E = Xmat(:,4);
    xE =  mT; xE(E == 1) = NaN; % turn GOOD index values into NaN so they don't print
    yE = dmT; yE(E == 1) = NaN; 
    
    if real_time == 0
        plot(yH,   'bx', 'MarkerSize',8); % Hausdorff
        plot(yE,   'g*', 'MarkerSize',8); % RJE    
        xlabel('Ordinal steps')
    elseif real_time == 1
        plot(xH,yH,   'bx', 'MarkerSize',8);
        plot(xE,yE,   'g*', 'MarkerSize',8);
        xlabel('Cumulative time (s) or displacement (cm)')
    end
    
    ylabel('Inter-STEP Int.')
    hold off
    
    % first order LOCATION transformation  
    subplot(nplots,1,3)
    ind = 1:numel(loc_trans);
    
    if real_time == 0
        plot(loc_trans,'k','Marker','.')
        hold on
        plot(ind(loc_trans > outlier_thr),loc_trans(loc_trans > outlier_thr),'r*')
        hold off
        xlabel('Successive events')
        
    elseif real_time == 1
        plot(mT, loc_trans, 'k','Marker','.') % use beg_ind:end_ind to account for ex_cap
        xlabel('Cumulative time (s) or displacement (cm)')
    end

    % y-axis label
    if strcmp(outlier_meth,'q')
        yalbel('Location quotiont')
    elseif strcmp(outlier_meth,'pc')
        ylabel('Location percent change')
    elseif strcmp(outlier_meth,'pd')
        ylabel('Location percent diff.')
    end
    
    % first order successive transformation  
    if do_suc == 1
        subplot(nplots,1,4)
        ind = 1:numel(suc_trans_all);

        if real_time == 0
            %plot(xind(beg_ind:end_ind), suc_trans_good(beg_ind:end_ind),'k','Marker','.')
            plot(suc_trans_all,'k','Marker','.')
            hold on
            plot(ind(suc_trans_all > outlier_thr), suc_trans_all(suc_trans_all > outlier_thr),'r*')
            hold off
            xlabel('Successive events')

        elseif real_time == 1
            plot(mT, suc_trans_all, 'k','Marker','.') % use beg_ind:end_ind to account for ex_cap
            xlabel('Cumulative time (s)')
        end

        % y-axis label
        if strcmp(outlier_meth,'q')
            yalbel('Successive quotiont')
        elseif strcmp(outlier_meth,'pc')
            ylabel('Successive percent change')
        elseif strcmp(outlier_meth,'pd')
            ylabel('Successive percent diff.')
        end
    else
        subplot(nplots,1,4)
        plot(nan)
    end
    %% manual exclusions 

    % Now we need to VISUALLY exclude any portion of the data
    % (independent of the other scrubbing methods)

    % interactively select the desired range(s)
    %
    %
    % ************
    % Note: although outliers are VISUALIZED, they are only EXCLUDED
    % after the following step is complete; skipping this step will
    % only exclude the end caps

    if ex_method == 2 || ex_method == 4

            % Use the cursor to mark the x-axis location(s), and press
            % RETURN when finished.

            man_ex = ginput; % x-values; will collect unlimited until ENTER is pressed

            man_ex = round(man_ex);
            man_ex = man_ex(:);
            
    end
         
end % make_plots

%% now concatenate the indices
inds = 1:size(Xmat,1);

H_all = Xmat(:,3);
E_all = Xmat(:,4);

% just the list of indices to exclude
H_ex = inds(H_all == 0); H_ex = H_ex(:);
E_ex = inds(E_all == 0); E_ex = E_ex(:);

% now we merge this
if ex_method == 0
   ex = ex_inds_import; % these will be indices, not an array of 1s and 0s; note: this will *not* include the end cap, but that is what we expect
elseif ex_method == 1
   ex = H_ex;
elseif ex_method == 2
   ex = [H_ex; man_ex];
elseif ex_method == 3
   ex = E_ex;
elseif ex_method == 4
   ex = [E_ex; man_ex];
elseif ex_method == 5
   ex = [];
end  

% now we delete extra points, if appropriate (the end cap exclusion has
% already been performed)

if del_itis == 0
    % don't adjust this
elseif del_itis == 1
    ex = [ex-1; ex; ex+1]; 
elseif del_itis == 2
    ex = [ex-2; ex-1; ex; ex+1; ex+2]; 
elseif del_itis == 3
    ex = [ex-3; ex-2; ex-1; ex; ex+1; ex+2; ex+3]; 
elseif del_itis >= 4
    ex = [ex-4; ex-3; ex-2; ex-1; ex; ex+1; ex+2; ex+3; ex+4]; 
end

ex = sort(unique(ex(ex>0))); % these are the indices to exclude (excludes "negative" indices)

% also need to delete numbers after the end of the file
ex = ex(ex<=numel(inds));

% a row for exporting
ex = ex'; 

% now get the structure for long_run_ts
final_ex = ones(size(H_all));
final_ex(ex) = 0; % so now included = 1 and excluded = 0

% and now this becomes a new, final column of Xmat
Xmat(:,5) = final_ex;
    
% 2013.10.12: having excluded all the outliers, let's just do a single pass
% run to calculate slope, preserving the "gaps" as  *empty space*

% xv = mT;  xv = xv(final_ex == 1);
% yv = dmT; yv = yv(final_ex == 1);
% 
% figure(300)
% plot(xv,yv)



%% after all exclusions (auto + manual), then we can find long runs
% we can look for runs that have a minimum number of elements

% 2013-09-12: we don't care about the merge_threshold right now, since all runs will
% by definition be close to the location. Besides, we calculate the global
% slope later on to identify linear non-stationarity in the data.

run_merge_meth = 'free'; % we don't use this param anymore since all Runs must be close to location 
run_merge_thr = Inf; % effectively ignores this parameter (for BEATS, but can re-instantiate it for Physionet)
%Xmat_input_to_long_run = Xmat

long_run_out = long_run_ts(Xmat, run_merge_meth, run_merge_thr, string_type, run_dur_thr, gap_dur_thr, do_slope); % will identify L and R alternations within this, after we select which is the "base" measure to look at

% now we choose whether we concatenate runs, or use the longest run

if strcmp(long_run_method,'s')
    % single (longest run)
    runs_ind = long_run_out.long_run_ind;
elseif strcmp(long_run_method,'c')
    % concatenate runs
    runs_ind = long_run_out.all_runs_ind;
end

% add this as a *new* column, replacing NaN with zero
runs_ind(isnan(runs_ind)) = 0;
Xmat(:,6) = runs_ind;

%% back to plot 

if plot_it == 1
    figure(50)
    subplot(nplots,1,2)
    hold on
    
    xO =  mT; xO(runs_ind == 1) = NaN; % exclude the *good* indices
    yO = dmT; yO(runs_ind == 1) = NaN;  

    if real_time == 0 
        plot(yO,   'ro', 'MarkerSize',8); % other cases
    elseif real_time == 1
        plot(xO,yO,'ro', 'MarkerSize',8); % other cases
    end

    hold off
end

%% some stats
stats.num_vals = numel(sort_ts);
stats.min_IEI = min(dmT);
stats.max_IEI = max(dmT);
stats.location = loc;
stats.H_inc_prc = 100 * sum(H_val) / numel(sort_ts); % this is just the method itself *not* other exclusions based on alternating steps or longest run!
stats.E_inc_prc = 100 * sum(E_val) / numel(sort_ts);
stats.loc_transf_max = loc_transf_max;
stats.suc_transf_max = suc_transf_max;
stats.runs_indices = runs_ind; % 1 = good, 0 = bad
stats.good_indices = find(runs_ind == 1); % the actual index lcoations
stats.ex_inds = ex;
stats.gap_dur_sum = long_run_out.gap_dur_sum;
stats.slope_pc_max_abs = long_run_out.slope_pc_max_abs;
stats.plot_slope_points = long_run_out.plot_slope_points;

    