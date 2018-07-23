function [output] = gait_sim(DTYPE,DPARAM,M,CV,SMTYPE,SMPARAM,N,HZ,ITER,PRC,PLOT_IT)

%
% function [output] = gait_sim(DTYPE,DPARAM,M,CV,SMTYPE,SMPARAM,N,HZ,ITER,PRC)
%
% create a noisy series, and see how outcome measures are affected by
% changes in "sampling rate"
%
% ** this version treats manipulations as REPEATED MEASURES as much as
% possible. For an independent variable version of this, see "gait_sim_bw"
%
% Order of operations
%   1. Run smooth_series ITER times, using max(N) values
%   2. Within each unique sequence, do the SMPARAM loop
%   3. Go from the CDF to the target values
%   4. From the "center" portion of the data, now do the N loop: the first N events, in sucession
%   5. After getting the desired N(n), then we recenter to target M 
%   6. Within each M loop, then we do the CV loop
%   7. Within each CV loop, then we do the HZ loop
%
%
% Figure has three plots:
%   1. Blue is original ITI series; red is the series with noise
%   2. Correlation of Original ITIs with Noisy ITIs
%   3. Poincare plot of the Original series

% calls function from:
% http://www.mathworks.com/matlabcentral/fileexchange/5091-generate-spatial-data/content/spatialPattern.m
%


%% hello

loc = which('gait_sim'); % note: must have a unique name and NOT an internal variable name
                      % in order to have the function called "beats", the
                      % internal variable must have a different name; i.e.,
                      % "beatsvar"
file_info = dir(loc);
save_date = file_info.date;

fprintf(['\n\n || Monte Carlo simulations of step time series, with outcome statistics and plots \n || Version: ' save_date '\n || For more information, see: \n ||   http://m3r.comp.nus.edu.sg/wiki/index.php/Main_Page \n ||   http://robjellis.net \n\n']);

%% initialize params

% using default values

if nargin >=8
    if M == 999;            
        %M =  [.40 .60 .90 1.35]; 
         M = [.125 .250 .500 1.00]; 
    end 
    
    if CV == 999;           
        CV = [1.25 2.5 5 10]; 
    end
    
    if SMPARAM == 999;   
        if strcmp(SMTYPE,'ma')
            SMPARAM = [3 5 9 17];
        elseif strcmp(SMTYPE,'wp')   
            %SMPARAM = [0 0.5 0.75 0.875 1.0]; 
            SMPARAM = [0 0.5 1.0];
        end
    end
    
    if N == 999;            
        N  = [25 50 100 200]; 
    end
    
    if HZ == 999;           
        HZ = [25 50 100 200 400 800 1600]; 
    end
    
    if PRC == 999;  
        PRC = 50:10:90; 
    end
else
    % we have a problem
    fprintf(' Error: not enough input parameters.\n\n')
end

if nargin < 8
   PLOT_IT = 0;
end


    % **********************************************
    % a few parameters (leave these here)
    do_dfa = 0; % will only work if series is long enough
    prc = PRC;
    %iter_prc = 95; % to check the 95th percentile of each round of iterations (assumes that all outcome stats are >= 0)
    loc_meth = 'mean';
    outlier_meth = 'pc';
    
    % **********************************************
    
    % how many rows do we have for outcome measures?
    ycheck = randn(100,1)+500;
    xcheck = cumsum(ycheck);
    num_events = Inf; % use the full series

    [td_check td_check_rows] = td_calcs_concat(xcheck, prc ,loc_meth, outlier_meth, do_dfa, num_events, 0); % for this iteration, these are the ground truth values 
    num_outcomes = size(td_check_rows,2);
    
%     % now we *add* a few extra cases, to track (1) min(Ygt), (2) max(Ygt), (3) max(ratio)
%     num_extra = 3;
%     num_outcomes = num_outcomes + num_extra;

%% outcome measures

% the FINAL structure will be set up to directly import into Statistica:
% the b/w subjects variables are columns; outcome measures are rows
% ("within-subjects")

%num_between = 5; 
%num_final_rows = numel(M) * numel(CV) * numel(SMPARAM) * numel(N) * numel(HZ) + 1; % +1 is becaus we have a row header
%num_final_cols = num_between + num_outcomes; % +num_between is because we have column headers

%td_mat_final = nan(numel(M),numel(CV),numel(SMPARAM),numel(N),numel(HZ),num_outcomes); % matrix format for other applications
%td_row_sign_final = nan(num_final_rows,num_final_cols);
%td_row_abs_final  = nan(num_final_rows,num_final_cols);

% 2014.03.06 - for the iRACE paper
% just save the 2.5, 50, and 97.5 percentile
irace_matrix = nan(3, numel(M), numel(CV), numel(N), numel(SMPARAM), num_outcomes); % percentiles, indiv plot x-axis, indiv plot sep lines, subplot rows, subplot cols, number of subplots

cur_row = 1; % cumulative counter

%% LOOP STARTS HERE

% define the integer seed values (from 1 to 10M)
SEED = randint(ITER,1,[1 10000000]);

if ITER > 1
    progressbar(0,0,0,0,0,0) 
end

for m = 1:numel(M)
    if ITER > 1
        progressbar((m-.0001)/numel(M),[],[],[],[],[])
    end
    for v = 1:numel(CV)
        if ITER > 1
        progressbar([],(v-.0001)/numel(CV),[],[],[],[])
        end
        for s = 1:numel(SMPARAM)
            if ITER > 1
            progressbar([],[],(s-.0001)/numel(SMPARAM),[],[],[])
            end
            for n = 1:numel(N)
                if ITER > 1
                progressbar([],[],[],(n-.0001)/numel(N),[],[])
                end

                for h = 1:numel(HZ)
                    if ITER > 1
                    progressbar([],[],[],[],(h-.0001)/numel(HZ),[])
                    end
                    
                    % get the holding bin for the ITER x outcome_measures
                    td_err_sign_tmp = nan(ITER,num_outcomes); % we want the first two variables to be rows and columnns
                    td_err_abs_tmp  = nan(ITER,num_outcomes); % we want the first two variables to be rows and columnns

                    %sim_vals_tmp = nan(ITER,num_extra); % for min(IEI), max(IEI), and max(ratio)
           
                    for i = 1:ITER
                        if ITER > 1
                        progressbar([],[],[],[],[],(i-.0001)/ITER);
                        end
                        
                        % get the ground truth (gt), smooth series, which
                        % will have EXACTLY the target mean and STD
                        % revised 2014.03.07 - we ALWAYS call a specific
                        % seed, so we can control the "repeated measures"
                        % nature of this process without having to
                        % re-engineer the entire loop structure
                        res = smooth_series(DTYPE, DPARAM, SEED(i), M(m), CV(v),  SMTYPE, SMPARAM(s), max(N), N(n), PLOT_IT); % we only look at gamma and weighted priors smoothings
                        
                        Ygt = res.val_smth_Nout; % these are the inter-event intervals; now need to get timestamps
                        Xgt = cumsum(Ygt);
                        
 
                        %% skip outlier stuff for now
                        do_outlier = 0;
                        
                        if do_outlier == 1
                            % let's run the outlier detector to see how the methods perform in terms of how much data they cut out
                            ex_cap = 4;
                            ex_method = 3; % use RJE auto method
                            ex_inds_import = [];
                            del_itis = 0;
                            long_run_method = 'c'; % perform concatenation
                            outlier_meth = 'pc';
                            transform_meth = 'pc';
                            outlier_thr = 50; 
                            run_dur_thr = 10; % in seconds
                            gap_dur_thr = 5; % in seconds
                            plot_out = 0;
                            do_slope = 0;
                            fnum = i;
                            fname = 'Simulation';

                            [outlier_Xmat outlier_stats] = td_outlier_detect(Xgt,[],ex_cap,ex_method,ex_inds_import,del_itis, long_run_method, loc_meth, outlier_meth, outlier_thr, run_dur_thr, gap_dur_thr, do_slope, plot_out,fnum,fname);
                        
                        
                            % here we get infor about how much of the data we retain from the Hausdorff vs Ellis methods, and
                            % also the largest observed ratio value >= median index

                            %sim_vals_tmp(i,1) = outlier_stats.min_IEI;
                            %sim_vals_tmp(i,2) = outlier_stats.max_IEI;
                            %sim_vals_tmp(i,3) = outlier_stats.E_max_rat;

                       end 
                        
                        %% Now we get the "ground truth" stats
                        % * td_calcs_cum can work with either one or two channels (e.g., L and R) of data
                        % * in the run below, Xgt is just a single column of data, so we don't worry about (1) channels or (2) outliers
                        
                            if ITER  == 1
                                make_td_plots = 1; % figures 150 and 155
                            else
                                make_td_plots = 0;
                            end
                         fignum = 150;   
                         use_interp = 0;
                         transform_meth = 'pc';
                         % do_dfa = 0 by default
                         
                        [td_gt td_gt_rows] = td_calcs_concat(Xgt, prc , loc_meth, transform_meth, do_dfa, num_events, make_td_plots,use_interp,fignum); % for this iteration, these are the ground truth values 

                        %dfa_val = td_gt.dfa_alpha;
                        
                        % also get the robust slope coordinates
                        %gt_slope_x = td_gt.xvals;
                        %gt_slope_y = td_gt.yvals;
                        
                        %% Now we get the (1) shorter length series and the (2) the downsampled series (ds)
                        % first, the shorter length
                        Xds = []; % just clear it to be safe
                        Yds = [];
                        
                        Xorig = Xgt(1:N(n)); % only need the x-values
                        
                        % now we down-sample that
                        sr = (1 / HZ(h)); % samples every xx *seconds*

                           Xds = ceil(Xorig / sr) * sr; % the new time stamp is rounded *up* to the next sampling point

                           Yds(1) = Xds(1);
                           Yds(2:numel(Xds)) = diff(Xds); % this is now the resampled ITI series, which is compared to original X series
                           
                           % get into a single column
                           Xds = Xds(:);
                           Yds = Yds(:);
                                                                                 
                            %% use td_calcs to do all the time-domain calculations 
                            % these are all performed, for each sequence, by td_calcs (which outputs in rows)

                            fignum = 155;
                            use_interp = 0;
                            [td_ds td_ds_rows] = td_calcs_concat(Xds, prc, loc_meth, transform_meth, do_dfa, num_events, make_td_plots,use_interp,fignum); % will be in a row


                            %% now we *compare* the down-sampled verion with the ground
                            % truth version using percent differences

                            % do we do this as SIGNED or ABSOLUTE ERROR? 
                            
                            % however, the *linear* difference may be more
                            % important, because (1) we always do ANOVA stats on PD vs HC, and (2) taking a percent
                            % difference of scores that are ALREADY percents can get even more confusing.
                            % Furthermore, RJE confirmed that this makes a huge difference in the "apparent" error.
                            % Finally, doing the simple signed or absolute difference is more directly relateable to the
                            % y-axis of a Bland-Altman plot
    
                            sign_meth = 'a'; % as of 2014.03.06, 'a' is still a good choice
                            
                            if strcmp(sign_meth,'a')
                                td_err_sign = (td_ds_rows - td_gt_rows); % sign is easier to interpret: "greater variability than ground truth" or "less variability than ground truth"
                            elseif strcmp(sign_meth,'p')
                                td_err_sign =  100 * (td_ds_rows - td_gt_rows) ./ td_gt_rows;
                            end
                            
                            % absolute error - not as useful!
                            td_err_abs = abs(td_err_sign);
                            
                            % put it into the temporary row
                            td_err_sign_tmp(i,:) = td_err_sign;
                            td_err_abs_tmp(i,:)  = td_err_abs;
                            
                            if PLOT_IT == 1
                                
                                interp_meth = 'v';
                                % interpolate the staircase of the DS series
                                isi_interp = staircase_interp(Yds,interp_meth);

                                isi_diff_interp = staircase_interp(diff(Yds),interp_meth);
                                
                                
                           plot_400 = 0;                          
                           if plot_400 == 1
                               % do same thing for abs difference
                               Ygt_dl = abs(Ygt - mean(Ygt));
                               Ygt_sd = abs(diff(Ygt));
                               Yds_dl = abs(Yds - mean(Yds));
                               Yds_sd = abs(diff(Yds));

                               figure(400)
                               subplot(1,2,1)
                               plot(sort(Ygt_dl),'b')
                               hold on
                               plot(sort(Yds_dl),'r')
                               hold off

                               subplot(1,2,2)
                               plot(sort(Ygt_sd),'b')
                               hold on
                               plot(sort(Yds_sd),'r')
                               hold off
                           end
                                
                          plot_ecdf = 0;
                          if plot_ecdf == 1
                              figure(175)
                              [fgt xgt] = ecdf(Ygt);
                              [fds xds] = ecdf(Yds);

                              plot(xgt,fgt,'b')
                              hold on
                              plot(xds,fds,'r')
                              hold off
                          end
%                                 % code by Zekun; RJE does this with KDE instead
%                                 numbins = round(numel(Yds)/10);
%                                 figure(175)
%                                 hist(Ygt,numbins);
%                                 hold on
%                                 if DTYPE == 'g'  
%                                     pd = fitdist(Yds,'Gamma');
%                                 elseif DTYPE == 'n' 
%                                     pd = fitdist(Yds,'Normal');
%                                 end
%                                 title('(From gait_sim.m)','Interpreter','none')
%                                 xlabel('Inter-event Intervals')
%                                 ylabel('Count')
%                                 
%                                 [bincounts,binpositions] = hist(Ygt,numbins);
%                                 binwidth = binpositions(2) - binpositions(1);
%                                 histarea = binwidth*sum(bincounts);
%                                 x_values = binpositions(1):0.001:binpositions(end);
%                                 fit = pdf(pd,x_values);
%                                 plot(x_values,histarea*fit,'m');
%                                 hold off
% 
%                                 Yfit = random(pd,size(Ygt));

                                % **********************
                                if PLOT_IT == 1
                                    extra = .01;
                                    xmin = min(Xds)*(1-extra);
                                    xmax =  max(Xds)*(1+extra);
                                    ymin = min(Yds)*(1-extra);
                                    ymax =  max(Yds)*(1+extra); 
                                    
                                    % make figures for easier use
                                    
                                    %% figure 101
                                    figure(101)
                                    plot(Xgt,ones(size(Xgt)),'bo')
                                    hold on
                                    plot(Xds,ones(size(Xgt)),'ro')
                                    hold off
                                    
                                    %% figure 102
                                    figure(102)
                                    subplot(2,1,1)
                                    plot(Xgt,Ygt,'b','LineWidth',2) % plot this behind
                                    hold on
                                    plot(Xds,Yds,'r','LineWidth',1) % downsampled (plot this in front)
                                    xlabel('Time (s)')
                                    ylabel('Inter-event Interval (s)')
                                    title('(From gait_sim.m)','Interpreter','none')

                                    axis([xmin xmax ymin ymax])
                                    hold off  
                                    
                                    subplot(2,1,2)
                                    % histogram
                                    hist(Ygt,30)
                                    xlabel('Inter-event interval (s)')
                                    ylabel('Count')
                                    
                                    % GT vs DS
                                    plot(Ygt,Yds,'.')
                                    xlabel('Original ITIs'); ylabel('Downsampled ITIs');
                                    
                                    %% figure 200        
                                    figure(200)

                                    subplot(1,11,1:3)
                                    plot(Xgt,Ygt,'b','LineWidth',2) % plot this behind
                                    hold on
                                    plot(Xds,Yds,'r','LineWidth',1) % downsampled (plot this in front)
                                    xlabel('Time (s)')
                                    ylabel('Inter-step Interval (s)')
                                    title('(From gait_sim.m)','Interpreter','none')

                                    axis([xmin xmax ymin ymax])
                                    hold off
                        
                                    %plot([ds_slope_x(1) ds_slope_x(2)],[ds_slope_y(1) ds_slope_y(2)],'Color','r','LineWidth',2); 
                                    %plot([gt_slope_x(1) gt_slope_x(2)],  [gt_slope_y(1) gt_slope_y(2)],'Color','b','LineWidth',1); 

                                    %plot([Xds(1) Xds(end)],[median(Yds) median(Yds)],'Color','m','LineWidth',1); % median
                                    %plot([Xgt(1) Xgt(end)],[median(Ygt) median(Ygt)],'Color','c','LineWidth',1); % median

                                    subplot(1,11,4:6)
                                    % sorted to make it easier to see
                                    % ** the sorted plot clearly shows why the percentile-based statistics are NOT ideal 
                                    % at low sampling rates ... mean-based methods are actually better

                                    plot(sort(Ygt),'b','LineWidth',2) % original: plot this behind
                                    hold on
                                    plot(sort(Yds),'r','LineWidth',1) % downsampled (plot this in front)
                                    %plot(sort(Yfit),'m','LineWidth',1); % from fitdist        
                                    plot(isi_interp.xi,    isi_interp.yiH,   'g','LineWidth',1) % hermite
                                    %plot(isi_interp.xi,
                                    %isi_interp.yiS,   'g') % cubic spline - doesn't work!

                                    % show the location of the interpolating points
                                    plot(isi_interp.x_disc,isi_interp.y_disc,'og')
                                    axis([-10 numel(Ygt)+10 ymin ymax]) 
                                    hold off
                                    xlabel('ISIs, sorted (s)')

                                    subplot(1,11,[7 8])
                                    for_box = [Ygt Yds];
                                    boxplot_rje(for_box,[2.5 97.5],0,0);
                                    axis([0.5 2.5 ymin ymax])
                                    xlabel('Box plots')
                                    hold off

                                    subplot(1,11,[9:11])
                                    plot(sort(diff(Ygt)),'b','LineWidth',2) % successive differences (ground truth series)
                                    hold on
                                    plot(sort(diff(Yds)),'r','LineWidth',1) % successive differences (downsampled series)
                                    plot(isi_diff_interp.x_disc,isi_diff_interp.y_disc,'og')
                                    plot(isi_diff_interp.xi,    isi_diff_interp.yiH,   'g','LineWidth',1) % hermite
                                    %plot(isi_diff_interp.xi,    isi_diff_interp.yiS,   'g') % cubic spline
                                    hold off
                                    %axis([-10 numel(Ygt)+10 ymin ymax])
                                    xlabel('diff(ISI), sorted (s)')
                                
                                end % plot figure 200

                                % RJE successive ratio illustration
                                plot_ratio = 0;
                                if plot_ratio == 1
                                    figure(225)
                                    rat_gt = Ygt(2:end) ./ Ygt(1:end-1);
                                    rat_ds = Yds(2:end) ./ Yds(1:end-1);
                                    plot(rat_gt,'b','LineWidth',2)
                                    hold on
                                    plot(rat_ds,'r','LineWidth',1)
                                    hold off
                                    xlabel('Successive IEI ratios')
                                    ylabel('Ratio value')
                                    title('(From gait_sim.m)','Interpreter','none')
                                end
                               
                               % uses figures(500) 
                               plot_kde = 0;
                               prctile_kde(Yds,'m',75,plot_kde)

                            end % PLOT_IT

                    end  % ITER loop
                    
                    %% now take percentiles and put into irace_matrix
                    prcs = prctile_nist(td_err_sign_tmp,[2.5 50 97.5]);
                    
                    for o = 1:num_outcomes
                        irace_matrix(:,m,v,n,s,o) = prcs(:,o);
                    end
                    
%                     %% now we get the Pth percentile of the iteration loop (assumes absolute error is taken).
%                     % note: prctile_nist assumes we have a 2D structure (N rows x M columns). 
%                     
%                     iter_prc = 95; % default, conservative
%                     
%                     tmp_prc = prctile_nist(td_err_abs_tmp,iter_prc); % yields a [1 x M column] vector for all outcome measures
%                     
%                     % for the observed distribution values, we have to do this a little differently
%                     vals_prc(1) = prctile_nist(sim_vals_tmp(:,1),5); % for min val
%                     vals_prc(2) = prctile_nist(sim_vals_tmp(:,2),95); % for max val
%                     vals_prc(3) = prctile_nist(sim_vals_tmp(:,3),iter_prc); % for ratio
                    
                    % we can visualize signed error
                    if ITER > 10 %&& PLOT_IT == 1
                        outlier_plot = 0;
                        figure(300)
                        clf % just to make sure we don't look at old results
                        boxplot_rje(td_err_sign_tmp,[2.5 97.5],300,outlier_plot); % use ";" to not print stats
                        title(['CV = ' num2str(CV(v)) '; N = ' num2str(N(n))])
                        axis([0 num_outcomes+1 -2 10])
                    end
                    
                    % note! this may not yield a linear relationship
                    % between, e.g., percentile values in PADM and the
                    % resulting tmp_prc value, because we are looking at
                    % *error* between the ground-truth and down-sampled
                    % versions
                    
                    %% store the results
                    % header row
                    
%                     if cur_row == 1;
%                         td_row_final(1, :) = [nan(1,num_between) 1:1:num_outcomes];
%                         td_row_final(2, :) = [M(m) CV(v) SMPARAM(s) N(n) HZ(h) tmp_prc vals_prc]; % single row, with appropriate column labels
%                         cur_row = 3; % next one starts here
%                     else
%                         % all other rows
%                         td_row_final(cur_row, :) = [M(m) CV(v) SMPARAM(s) N(n) HZ(h) tmp_prc vals_prc]; % single row, with appropriate column labels
%                         cur_row = cur_row + 1; % increment
%                     end
                    
                    
                    % store in the multi-dimensionsal matrix                           
                    %td_mat_final(m,v,s,n,h,:) = [tmp_prc vals_prc];
                    
                end % HZ loop
           end % N loop     
        end % SM loop           
    end % v (CV) loop
end % m (MEAN) loop

if ITER >= 1
    progressbar(1,1,1,1,1,1)
end
%% plots

%boxplot_rje(td_row_final(:,6:end));

% plot for the iRACE paper
nn = numel(N);
ss = numel(SMPARAM);
if ITER > 100
    for o = 1:num_outcomes
       figure(400+o)
       ctr = 1;

       % figure out y axis limits
       all_o = irace_matrix(:,:,:,:,:,o);
       all_o = all_o(:); % for simplicity
       ymin = min(all_o);
       ymax = max(all_o);

       for s = 1:ss
           for n = 1:nn
                subplot(ss,nn,ctr)
                    for v = 1:numel(CV)
                        data = irace_matrix(:,:,v,n,s,o); % will be an N x M simple matrix
                            if v == 1
                               color = 'r';
                            elseif v == 2
                               color = 'g';
                            elseif v == 3
                               color = 'b';
                            elseif v == 4
                               color = 'k';
                            end

                            for p = 1:3 % percentiles (just the outer ones)
                                plot(data(1,:),color)
                                hold on
                                plot(data(3,:),color)
                                if numel(M) > 1
                                    axis([1 numel(M) ymin ymax])
                                else
                                    axis([0 2 ymin ymax])
                                end

                                if n == 1 % only for the left subplot of each row
                                    ylabel('Signed error')
                                elseif n > 1
                                    set(gca,'YTickLabel',{'','','',''}) % don't display label
                                end

                                if s < 3
                                    set(gca,'XTickLabel',{'','','',''}) % don't display label
                                elseif s == 3 % only for the bottom row
                                    xlabel('IEI Mean')
                                    set(gca,'XTickLabel',{'.125','.250','.500','1.00'})
                                end
                            end                       
                    end
                    hold off % after all CV lines are drawn
                ctr = ctr + 1;           
           end
       end
    end
end
%% outputs
    output.M =          M;
    output.CV =         CV;
    output.SMPARAM =    SMPARAM;
    output.N =          N;
    output.HZ =         HZ;
    output.ITER =       ITER;
    output.num_outcomes = num_outcomes;
    % just copy the last sequence made for easy access
    output.Ygt = Ygt;
    output.Yds = Yds;
    %output.td_calc_prc = prc;
    %output.iter_prc = iter_prc;
    
    output.irace_matrix = irace_matrix;
    
    %output.td_mat_final = td_mat_final;
    %output.td_row_final = td_row_final;
    
end
