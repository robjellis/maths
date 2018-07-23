function [output output_rows] = td_calcs_concat(Xmat, prc, loc_meth, transform_meth, do_dfa, max_events, plot_it, use_interp, fignum)
% 
% function output = td_calcs_concat(Xmat, prc, loc_meth, transform_meth,
% do_dfa, plot_it)
%
% a version of td_calcs that allows *multiple* runs to be concatenated
% first
% before performing time domain calculations
%
% Xmat values are *timestamps*, not inter-event intervals. It will either be:
% - 1 columns: BEATS data. We have already found a *single* continuous segment of data, and assumes we don't need to worry about outliers
% - 2 columns: HR data. col1 = timestamps, col2 = valid (1) vs invalid (0) timestamps
% - 3/4 columnns: Gait data. co1 = timestamps; col2 = channel (L = 1 and R = 2); col3 = Hausdorff method; col4 = rje outlier identification (0 = outlier)
%
% "loc_meth" is either 'l' for lambda (calls cent_tend.m by RJE), 'kde'
% (external function kde.m),  or 'mh' for midhinge. Also, if an actual value is input (e.g., .542), we will use that directly as loc.
%
% "transform_meth" is either 'ad' for absolute diffrence, 'ld' for log difference, 'pc' for percent change, or 'pd' for absolute percentage difference
% 
% "prc" is a vector of desired percentiles; e.g., [50 60 70]

% "max_events" allows a threshold to only use the first N events (used for iRACE and iMOTIV calculations)
%
% measures calculated:
% - ** note: PCI is calculated in a separate function (pci_calc.m)
% - SD and CVSD
% - RMS and CVRMS: root mean square successive differences
% - PCI in the original formulation and a modified one by RJE
% - PQL: percentile of quotients relative to the location
% - PSQ: percentile of successive quotients
% - DFA alpha (if do_dfa == 1): per standard Hausdorff formulations; will return as NaN if the series isn't long enough
%
% rje | version = 2013.09.13
%

if nargin < 2
    %prc = input('\n Enter the desired percentiles; e.g., [10 50 90]: ');
     prc = [10 50 90];
end

if nargin < 3
    % assume the simplest case: uses mean
    loc_meth = 'mean';
end

if nargin < 4
    transform_meth = 'ad';
end

if nargin < 5
   do_dfa = 1;
end

if nargin < 6
    max_events = Inf;
end

if nargin < 7
    plot_it = 0;
end

if nargin < 8
    use_interp = 0;
end

if nargin < 9
    fignum = 150;
end

%% intro

prc = prc(:);

%% get the locations of the run(s) (always from the last column)

% note: if we input a *single* vector of numbers, it is assumed that ALL
% data points are valid (i.e., we just extract the longest stable segment
% for the BEATS project). So, we need a little help here ...

if min(size(Xmat)) == 1;
   run_inds = ones(size(Xmat,1),1);
else
   run_inds = Xmat(:,end);
end

run_inds = [0; run_inds; 0]; % the zeros allow us to identify beginning ([0 1 1 1]) or ending (e.g., [1 1 1 0]) runs
diffrun = diff(run_inds); % indexes will be correct to mark timestamps
    
startIndex = find(diffrun == 1); % find ALL changes from 0 to 1

endIndex = find(diffrun == -1); % final ALL changes from 1 to 0, or a 1 at end of vector
        
endIndex = endIndex - 1; % the last 1 before a 0; the "-1" makes this work correctly since we are using diff()

nruns = numel(startIndex); % so we can loop through this

%% output variable set up
if nruns == 0
    % no valid data
    
    output.num_total =      NaN;
    output.num_valid_orig = NaN;
    output.num_valid_capped = NaN;
    output.slope_pc_concat = NaN;
    output.num_incl =       NaN;
    output.duration =       NaN;
    output.percent_valid =  NaN;
    output.mean =           NaN;
    output.median =         NaN;
    output.trimean =        NaN;
    output.sd    =          NaN;
    output.cvsd  =          NaN;
    output.ad90 =           NaN; % 90th percentile of absolute deviations from the median
    output.cvad90 =         NaN;
    output.rms   =          NaN;
    output.cvrms =          NaN;
    output.slope_drift =    NaN; 
    output.xvals =          NaN;
    output.yvals =          NaN;
    output.dfa_alpha =      NaN;
    output.percentiles =    NaN;
    
    output.loc_transf_mn =  NaN;
    output.loc_transf_max = NaN;
    output.loc_prc =        NaN;
    output.lm_mean_loc =    NaN;
    output.um_mean_loc =    NaN;
    output.dm_mean_loc =    NaN;
    
    output.suc_transf_mn =  NaN;
    output.suc_transf_max = NaN;    
    output.suc_prc =        NaN;
    output.lm_mean_suc =    NaN;
    output.um_mean_suc =    NaN;
    output.dm_mean_suc =    NaN;
    
    output.phiR_num = NaN;
    output.phiL_num = NaN;
    output.phiM_num = NaN;
    
    output.phiR_abs_prc = NaN;
    output.phiL_abs_prc = NaN;
    output.phiM_abs_prc = NaN;
    
    output.phiR_cv = NaN;
    output.phiL_cv = NaN;
    output.phiM_cv = NaN;
    
    output.pciR = NaN;
    output.pciL = NaN;
    output.pciM = NaN;
    
    output_rows = NaN(1,7); % as of 2014.03.06
    
    return % skip the whole rest of the process to save time
end

%% a few calculations on the full segment

Xfull = Xmat(:,1); 
Yfull = diff(Xfull);

Yfull = [NaN; Yfull];
Xfull(1) = NaN; % to match Yfull

Yvalid = Xmat(:,end); % the last column will always be the *FINAL* set of included (= 1) and excluded (= 0) cases

% now we just turn excluded indices to NaN; diff will also work fine
X = Xfull;
Y = Yfull;

Y(Yvalid == 0) = NaN;
X(Yvalid == 0) = NaN;

% make a copy of this, with NaNs preserved, for better plotting
Xplot = X;
Yplot = Y;

   
% this is only valid for doing (1) location and (2) slope; 

X(isnan(X))         = [];
Y(isnan(Y))         = [];
    
    

% **************************
%% PCI (Plotnik et al)

% for simplicity, this has all been moved to a separate script:
% pci_calc_concat.m, since it is not relavant to always do PCI (e.g., for
% heart rate data or tempo data)

% revised: we justify PCI for a single column of gait data (simulated) by
% taking alternating events (done within the script)
plot_phi = plot_it;
min_L_R = 3; % at least 3 L and R events in a string
toe_off = [];
phi_max_ev = 0.5*max_events - 1; % assuming we have N valid intervals, we get a max of 0.5*N - 1 phi values for L and R

pci_res = pci_calc_concat(Xmat, plot_phi, phi_max_ev, min_L_R, toe_off);
   
%% now the loop to concatenate runs!

% start with [] and concatenate
loc_transf_concat = [];
suc_transf_concat = [];
sq_vals_concat    = []; % we will do RMSSD on this after concatenation
Yrun_concat       = []; % this is needed for RMSSD

loc = NaN; % placeholder

for n = 1:nruns
    loc_trans = []; % reset
    suc_trans = []; % reset
    slope_pc_all = nan(nruns,1);
    
    % Xrun and Yrun *do* include gaps, so we can't do slope calculations here!
    Xrun = Xfull(startIndex(n):endIndex(n)); % timestamps
    Yrun = diff(Xrun);

        
    npoints = 2^10; % used by kde methods
    
    if ischar(loc_meth)
        if strcmp(loc_meth,'bf')         
        % brute force method: find the value from which the most number of values have a percentage deviation of <= thr

        data_prc = 0; % ignore this rule for now
        brute_out = cent_tend_brute(Yrun,outlier_thr,data_prc);

        loc = brute_out.location;
        elseif strcmp(loc_meth,'kde_matlab')

            % simple rule: if we don't have enough unique, then just use mode
            prc_unique = 100 * numel(unique(Yrun)) / numel(Yrun);

            if prc_unique < 50
                % just use mode to be safe
                loc = mode(Yrun);
            else

                range = max(Yrun) - min(Yrun);
                MIN = min(Yrun) - range/4; 
                MAX = max(Yrun) + range/4;

                [f,xi] = ksdensity(Yrun,linspace(MIN,MAX,npoints)); 

                loc = min(xi(f == max(f))); 
            end
        elseif strcmp(loc_meth,'kde_botev')

            [bandwidth,density,xmesh] = kde(Yrun,npoints);

            loc = min(xmesh(density == max(density))); % take min so we get the lowest peak

        elseif strcmp(loc_meth,'tm') 
            % http://en.wikipedia.org/wiki/Trimean

            % b. for Physionet, we use the trimean rather than the median, since
            % it will better reflect the positive skew of the distribution

            loc = (prctile_nist(Yrun,25) + 2* prctile_nist(Yrun,50) + prctile_nist(Yrun,75)) / 4; % Tukey formulation

        elseif strcmp(loc_meth,'med')
            % regular median
            Ytmp = Yrun(isnan(Yrun) == 0);
            loc = median(Ytmp);

        elseif strcmp(loc_meth,'mean')
            Ytmp = Yrun(isnan(Yrun) == 0);
            loc = mean(Ytmp);
        end
    else
        % if we use a number as input, then we use this number as the location
        loc = loc_meth;
    end
    
    %% slope calculations
    % this is done within long_run_ts now
    
    
    %% outcome measure transformations (within each run)
    % note: we don't do any outlier detection here; assume outliers have
    % already been detected as part of td_outlier_detect
    %
    % Gaps are *included* in Yrun!
    
    % 2013.10.11: td_outlier_detect now does a two-step (deviations from
    % location first, followed by successive deviations); here, we can do
    % them at the same time, since we only deal with RUNS
   
    %%%%%%%%%%%%%%%%
    % rmssd since it is easy and doesn't use loc, just mean from within each run
    
    sq_vals = diff(Yrun) .^ 2;
        
    %%%%%%%%%%%%%%%
    % the other calculations are a bit more complex
    
    if strcmp(transform_meth,'adl')
        % simple absolute difference
        loc_trans = abs(Yrun - loc); % note: this subtracts out the local mean
        suc_trans = abs(diff(Yrun));
        
    elseif strcmp(transform_meth,'adg') 
        % we use global estimate of location after we get Yrun
        loc_trans = NaN; % we use Yrun; this is just a place holder
        suc_trans = abs(diff(Yrun));        
    elseif strcmp(transform_meth,'ld')
        
    elseif strcmp(transform_meth,'pc') % percent change from previous beat / from location     
        % deviation from location
        loc_trans = 100 * abs(Yrun - loc) / loc;

        % successive 
        for i = 2:numel(Yrun)
            suc_trans(i-1) = 100 * abs(Yrun(i) - Yrun(i-1)) / Yrun(i-1);
        end 
        
    elseif strcmp(transform_meth,'pd') % percent difference  
        % deviation from location
        loc_trans = 100 * abs(Yrun - loc) / (0.5 * (Yrun + loc));


        % successive 
        for i = 2:numel(Yrun)
            suc_trans(i-1) = 100 * abs(Yrun(i) - Yrun(i-1)) / (0.5 * (Yrun(i) + Yrun(i-1)));
        end        
        
    elseif strcmp(transform_meth,'q')
        % deviation from location
        for i = 1:numel(Yrun)
            if Yrun(i) >= loc
                loc_trans(i) = Yrun(i) / loc;
            else
                loc_trans(i) = loc / Yrun(i);
            end
        end

        % successive
        for i = 2:numel(Yrun)
            if Yrun(i) >= Yrun(i-1)
                suc_trans(i-1) = Yrun(i) / Yrun(i-1);
            else
                suc_trans(i-1) = Yrun(i-1)/ Yrun(i);
            end
        end

    end 

    loc_trans = loc_trans(:);
    suc_trans = suc_trans(:);
    

    %% now do the actual concatenation
    
    loc_transf_concat = [loc_transf_concat; loc_trans];
    suc_transf_concat = [suc_transf_concat; suc_trans]; 
    sq_vals_concat   = [sq_vals_concat; sq_vals];
    Yrun_concat      = [Yrun_concat; Yrun];
    
  
    
end % nruns loops

if strcmp(transform_meth,'adg')
    Yrun_concat(isnan(Yrun_concat)) = []; 
    % for method 'adg', we need to now finally subtract global mean

    loc_transf_concat = abs(Yrun_concat - mean(Yrun_concat)); % can do trimean or something else later
end

% get the maximum abs. slope across runs
slope_pc_max = max(slope_pc_all);

%% after the concatenation loop is done, now do calculations
% note: we already terminate td_calcs_concat if there are no runs

% 0. general data clean-up

    % get rid of NaN
    loc_transf_concat(isnan(loc_transf_concat)) = [];
    suc_transf_concat(isnan(suc_transf_concat)) = [];
    sq_vals_concat(isnan(sq_vals_concat)) = [];
    Yrun_concat(isnan(Yrun_concat)) = []; 

    % save the original number of events
    Yrun_concat_orig = Yrun_concat;
             
    
%% 1a. Take the first N events
% subtract 1 so RMSSD etc. work 

if max_events == Inf
    % keep all events; don't need to redefine
    
else
    if max_events < numel(loc_transf_concat)
        loc_transf_concat = loc_transf_concat(1:max_events);
    %else
        %loc_transf_concat = NaN; % just keep the full series
    end
    
    if max_events < numel(suc_transf_concat)
        suc_transf_concat = suc_transf_concat(1:max_events);
    %else
        %suc_transf_concat = NaN;
    end
    
    if max_events < numel(sq_vals_concat)
        sq_vals_concat = sq_vals_concat(1:max_events);
    %else
        %sq_vals_concat = NaN;
    end
    
    if max_events < numel(Yrun_concat)
        Yrun_concat = Yrun_concat(1:max_events);    
    %else
        %Yrun_concat = NaN;
    end
end

    % now let's find the simple linear slope using polyfit
    % update: only do this for the CAPPED events!
    
    % need to have at least 10 events
    if numel(Yrun_concat) >= 10
        xv = 1:numel(Yrun_concat); xv = xv(:); 
        yv = Yrun_concat; yv = yv(:); 
        b = polyfit(xv,yv,1); % linear 

        y2 = b(1) * xv(end) + b(2);
        y1 = b(1) * xv(1)   + b(2);

        % percentage change relative to the mean
        slope_pc_concat = 100 * (y2 - y1) / mean(Yrun_concat);
    
        if plot_it == 1
            figure(200)
            plot(Yrun_concat_orig)
            hold on
            % now just plot the slope for the first n events
            plot([xv(1) xv(end)],[y1 y2],'r')
            hold off
        end
    else
        slope_pc_concat = NaN;
        plot_it = 0; % cancel the plot!
    end
    
    nloc = numel(loc_transf_concat); % will be "1" if max_events condition is not met
    nsuc = numel(suc_transf_concat);

% 1b. sort the data - valid for all measures
    sort_loc_transf_concat = sort(loc_transf_concat);
    sort_suc_transf_concat = sort(suc_transf_concat);


% interpolate
if use_interp == 1
    interp_meth = 'v';
    loc_interp = staircase_interp(loc_transf_concat,interp_meth);
    suc_interp = staircase_interp(suc_transf_concat,interp_meth);
end


% 2a. standard measures on Yrun_concat
sd    = std(Yrun_concat);
md    = median(Yrun_concat);
cvsd  = 100 * sd / mean(Yrun_concat); % note: if mean(Yrun_concat) approaches zero, then CV will be very large
rms   = sqrt(mean(sq_vals_concat));
cvrms = 100 * rms / mean(Yrun_concat); % rms relative to the mean (similar to coefficient of variation)

% 90th percentile of absolute deviations from the median
ad90 = prctile_nist(abs(Yrun_concat - md),90);
cvad90 = 100 * ad90 / md;

% mean and median absolute deviation (mad_rje is a separate function by RJE)
% update: don't worry about this now
xx = mad_rje(Yrun_concat);
mnad = xx.mnad;
cvmnad = xx.cvmnad;
mdad = xx.mdad;
cvmdad = xx.cvmdad;

% 3a. mean value of loc_prc and suc_prc: we *don't* do a "CV" version of this
loc_transf_mn = mean(loc_transf_concat);
suc_transf_mn = mean(suc_transf_concat);

% 3b. percentiles of the actual data
    loc_prc = prctile_nist(loc_transf_concat,prc);
    loc_prc = loc_prc'; % a row

    suc_prc = prctile_nist(suc_transf_concat,prc);
    suc_prc = suc_prc'; % a row

% 3c. max value of loc_transf_concat and suc_transf_concat

% ** done by td_outlier_detect

% 3d. "upper mean": mean of the upper P percent of all datapoints

% note for upper and lower mean: for local deviations, because we calculate the location for *each* run separately, 
% when we concatenate values, we won't see the "bins" in the same way that we will for absolute successive diffs.

    % get the index locations
    um_inds_loc = ceil(prc/100 * nloc + .10e-10);
    um_inds_suc = ceil(prc/100 * nsuc + .10e-10);
    
    % adjust if needed
    um_inds_loc(um_inds_loc < 1) = 1;
    um_inds_suc(um_inds_loc < 1) = 1;
    
    for p = 1:numel(prc)
        um_mean_loc(p) = mean(sort_loc_transf_concat(um_inds_loc(p):end));
        um_mean_suc(p) = mean(sort_suc_transf_concat(um_inds_suc(p):end));
    end
    
% 4. "lower mean": mean of lower P percent of all datapoints (uses same PRC
%    values as upper mean; we chose which of the two measures is more of
%    interest for us

    % get the index locations
    lm_inds_loc = floor(prc/100 * nloc - .10e-10);
    lm_inds_suc = floor(prc/100 * nsuc - .10e-10); 
    
    % adjust if needed
    lm_inds_loc(um_inds_loc < 1) = 1;
    lm_inds_suc(um_inds_loc < 1) = 1;

    for p = 1:numel(prc)
        lm_mean_loc(p) = mean(sort_loc_transf_concat(1:lm_inds_loc(p)));
        lm_mean_suc(p) = mean(sort_suc_transf_concat(1:lm_inds_suc(p)));
        
    end

        % now the difference between each pair
        dm_mean_loc = um_mean_loc - lm_mean_loc;
        dm_mean_suc = um_mean_suc - lm_mean_suc;
        
% **************************
% DFA alpha (requires external function, with modifications made by RJE
% http://www.mathworks.com/matlabcentral/fileexchange/19795-detrended-fluctuation-analysis

% this is also done on the cleaned up series for fairness

if do_dfa == 1 && numel(Yrun_concat) > 10 % second argument saves us from doing DFA on NaN
    
    dfa_meth = 3; % 1 = Peng, Hausdorff; 2 = Damouras; 3 = Jordan
    
    if dfa_meth == 1
        % 4 to 60 is standard for gait analyis (rje checked Hausdorff papers)
        [D alpha] = DFA_main(Yrun_concat,4,60,20,'g'); % Yrun_concat is the inter-strike interval series

    elseif dfa_meth == 2
        % below is based on Damouras; more accurate?
        [D alpha] = DFA_main(Yrun_concat,16,floor(numel(Yrun_concat)/9),[],'d'); 
    elseif dfa_meth == 3
        % based on Jordan 2006
        [D alpha] = DFA_main(Yrun_concat,4,floor(numel(Yrun_concat)/4),[],'d'); 
    end
    dfa_alpha = alpha; 
else
    dfa_alpha = NaN;
end


%% figures
    
    if plot_it == 1  
        figure(fignum)
        if numel(loc_transf_concat) > 1
            % then we do 4 rows
            r = 4;
        else
            % two rows
            r = 2;
        end
        
        subplot(r,2,1) % original
        plot(Xfull,Yfull,'Marker','.','Color',[0 0 0]) 
        xlabel('Cumulative time (s)')
        ylabel('Original IEI (s)')
        
        subplot(r,2,2) % get rid of outliers in plot
        plot(Xplot,Yplot,'Marker','.','Color',[0 0 0]) % NaNs will be gaps
        xlabel('Cumulative time (s)')
        ylabel('Outlier-scrubbed IEIs')
        
        %cv_calc = 100 * nanstd(Yplot) / nanmean(Yplot)
        
        % sort the actual IEIs to highlight bins
        
        subplot(r,2,3)
        plot(sort(Yfull),'Marker','.','Color',[0 0 0]) 
        xlabel('Sorted IEIs')
        ylabel('Original IEIs')
        
        subplot(r,2,4)
        plot(sort(Yplot),'Marker','.','Color',[0 0 0]) 
        xlabel('Sorted IEIs')
        ylabel('Outlier-scrubbed IEIs')
        
                % adjust the axes so we can *really* see the bins
                bin_size = diff(sort(Yplot));
                bin_size = bin_size(bin_size > .001); % ignore subtle bins
                bin_size = median(bin_size); % to be representative


                % how do we round?
                round_to = 1 / bin_size;

                % now determine the min and max values
                ymin = floor(min(sort(Yplot))*round_to) / round_to;
                ymax =  ceil(max(sort(Yplot))*round_to) / round_to;
                axis([0 numel(Yplot)+1 ymin ymax])
        
        % transformation
        
        if numel(loc_transf_concat) > 1
            subplot(r,2,5) % concatenated deviations
            plot(loc_transf_concat,'Marker','.','Color',[0 0 0])
            xlabel('Concatenated events')
            ylabel('Deviation transformation')

            subplot(r,2,6) % concatenated successive
            plot(suc_transf_concat,'Marker','.','Color',[0 0 0])
            xlabel('Concatenated events')
            ylabel('Successive transformation')

            % sorted IEIs to make the percentile measure clear

            % deviations
            subplot(r,2,7)
            plot(sort_loc_transf_concat,'Marker','.','Color',[0 0 0]) 
            hold on
            if use_interp == 1
                plot(loc_interp.xi,loc_interp.yiH,'b')
            end
            axis([0 numel(sort_loc_transf_concat)+1 0 max(sort_loc_transf_concat)])
            xlabel('Sorted concatenated events')
            ylabel('Sorted version of above')

            % plot the percentile locations
            prc_ind_loc = prc/100 * nloc; % this is exact

            for p = 1:numel(prc)
                plot([prc_ind_loc(p) prc_ind_loc(p)],[min(sort_loc_transf_concat) max(sort_loc_transf_concat)],'r')
            end
            hold off

            % successive
            subplot(4,2,8)
            plot(sort_suc_transf_concat,'Marker','.','Color',[0 0 0]) 
            hold on
            if use_interp == 1
                plot(suc_interp.xi,suc_interp.yiH,'b')
            end
            axis([0 numel(sort_suc_transf_concat)+1 0 max(sort_suc_transf_concat)])
            xlabel('Sorted concatenated events')
            ylabel('Sorted version of above')

            prc_ind_suc = prc/100 * nsuc; % this is exact

            for p = 1:numel(prc)
                plot([prc_ind_suc(p) prc_ind_suc(p)],[min(sort_suc_transf_concat) max(sort_suc_transf_concat)],'r')
            end

            hold off
        end
 
% This method won't help us ...        
%         % use KDE to get the CDF
%         subplot(4,2,7)
%         [bandwidth,density,xmesh,cdf]=kde(loc_transf_concat);
%         plot(xmesh,cdf,'k');
%         xlabel('IEI deviation transformation')
%         ylabel('KDE of PDF')
%         axis([min(xmesh) max(xmesh) 0 1])
%         
%         subplot(4,2,8)
%         [bandwidth,density,xmesh,cdf]=kde(suc_transf_concat);
%         plot(xmesh,cdf,'k');
%         xlabel('IEI successive transformation')
%         ylabel('KDE of PDF')        
%         axis([min(xmesh) max(xmesh) 0 1])
        
 
    end % plot_it

    
%% outputs
% note: output will be NaN if there is no valid data read in
%
% "Yrun_concat" is just the GOOD data (after we cap the number of events)
    output.num_total = numel(Xfull) - 1; % -1 to count inter-event intervals
    num_valid_orig = numel(Yrun_concat_orig(Yrun_concat_orig > 0)); % will be "0" if we only have NaN
    
    output.num_valid_orig = num_valid_orig;
    output.num_valid_capped = numel(Yrun_concat(Yrun_concat > 0));
    output.slope_pc_concat = slope_pc_concat;
        
    output.duration = sum(Yrun_concat);
    output.percent_valid = 100 * sum(Yrun_concat) / nansum(Yfull); % percent of the total recording
    output.mean = mean(Yrun_concat);
    output.median = median(Yrun_concat);
    trimean = (prctile_nist(Yrun_concat,25) + 2 * prctile_nist(Yrun_concat,50) + prctile_nist(Yrun_concat,75)) / 4;
    output.trimean = trimean;
    output.location = loc; % just provide the final location
    output.sd    = sd;
    output.cvsd  = cvsd;     
    output.ad90  = ad90;
    output.cvad90 = cvad90;
    output.rms   = rms;
    output.cvrms = cvrms;    
    %output.mnad = mnad;
    %output.cvmnad = cvmnad;
    %output.mdad = mdad;
    %output.cvmdad = cvmdad;
    output.slope_pc_max = slope_pc_max;
    %output.xvals = xvals;
    %output.yvals = yvals;
    output.dfa_alpha = dfa_alpha; 
    output.percentiles = prc;    % will be in rows
    output.loc_transf_mn = loc_transf_mn;
    output.loc_prc = loc_prc;
    output.lm_mean_loc = lm_mean_loc;
    output.um_mean_loc = um_mean_loc;
    output.dm_mean_loc = dm_mean_loc;
    output.suc_transf_mn = suc_transf_mn;
    output.suc_prc = suc_prc;
    output.lm_mean_suc = lm_mean_suc;
    output.um_mean_suc = um_mean_suc;
    output.dm_mean_suc = dm_mean_suc;
    
    output.phiR_all = pci_res.phiR;
    output.phiL_all = pci_res.phiL;
    output.phiM_all = pci_res.phiM;
    
    output.phiR_num = pci_res.phiR_num;
    output.phiL_num = pci_res.phiL_num;
    output.phiM_num = pci_res.phiM_num;
    
    output.phiR_mn = pci_res.phiR_mn; 
    output.phiL_mn = pci_res.phiL_mn;
    output.phiM_mn = pci_res.phiM_mn;
    
    output.phiR_abs_prc = pci_res.phiR_abs_prc;
    output.phiL_abs_prc = pci_res.phiL_abs_prc;
    output.phiM_abs_prc = pci_res.phiM_abs_prc;
    
    output.phiR_sd = pci_res.phiR_sd;
    output.phiL_sd = pci_res.phiL_sd;
    output.phiM_sd = pci_res.phiM_sd;
    
    output.phiR_cv = pci_res.phiR_cv;
    output.phiL_cv = pci_res.phiL_cv;
    output.phiM_cv = pci_res.phiM_cv;
    
    output.phiR_altcv = pci_res.phiR_altcv;
    output.phiL_altcv = pci_res.phiL_altcv;
    output.phiM_altcv = pci_res.phiM_altcv;    
    
    output.pciR = pci_res.pciR;
    output.pciL = pci_res.pciL;
    output.pciM = pci_res.pciMer;


    % for Physionet HRV
    % REVISED on 2014.03.06 for the iRACE paper (only a small set of stats output)
    %output_rows = [output.num_incl output.mean cvsd cvrms loc_transf_mn um_mean_loc suc_transf_mn um_mean_suc]; 
    output_rows = [output.mean cvsd cvrms pci_res.phiM_cv]; 
