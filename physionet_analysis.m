function [output] = physionet_analysis
% extract signal from physionet signals, based on the database
%
% 1: gaitdb (separate left and right timestamps)
% 2:
% 3: 
%
% updated = 2013.08.16

%% hello

loc = which('physionet_analysis'); % note: must have a unique name and NOT an internal variable name

file_info = dir(loc);
save_date = file_info.date;

fprintf(['\n || Signal processing / event extraction pipeline for Physionet datasets \n || Version: ' save_date '\n || For more information, see: http://robjellis.net \n\n']);


%% What type of dataset are we using?
db = input([' Select the database to analyze:'...
            '\n Physionet gait databases \n   [1] gaitpdb <ENTER> \n   [2] gaitdb \n   [3] gaitndd \n   [4] umwdb ',...
            '\n Pre-processed (strike detected) gait data \n   [11] Hausdorff method for gaitpdb',...          
            '\n Physionet R-R interval databases \n   [21] xxxx ',...
            '\n Physionet EKG databases  \n   [31] mvtdb ',...
            '\n Simple *.txt file inputs \n   [91] one channel, timestamps   \n   [92] one channel, inter-event intervals \n   [93] xxx'...
            '\n   --> ']);

if isempty(db)
    db = 1;
end

% Gait: http://physionet.org/physiobank/database/#gait
% IBI: http://physionet.org/physiobank/database/#rr

% select files, based on the database

if db == 1
   % only want the *_01.txt walks
   filt = '*_01.txt';
   grpE = 'Pt';
   grpC = 'Co';
elseif db == 2
   filt = '*.txt';    
elseif db == 3
   filt = '*.txt';
elseif db == 4
   filt = '*.txt';
elseif db == 11
   filt = '*.mat';
   grpE = 'Pt';
   grpC = 'Co';
elseif db == 21
   filt = '*.txt';   
elseif db >= 91
   filt = '*.txt';  
end

if db < 10
    [filenames, pathname] = uigetfile( {filt,'Physionet files'},'Select desired files','MultiSelect', 'on');
elseif db < 20
    [filenames, pathname] = uigetfile( {filt,'Physionet files'},'Select desired files','MultiSelect', 'off');
elseif db < 30
   % Physionet IBI files 
elseif db < 90
   % raw EKG files
elseif db < 100
    [filenames, pathname] = uigetfile( {filt,'Text files'},'Select desired files','MultiSelect', 'on');
end

% how many files?
if ischar(filenames)
   numf = 1;
   dt = 1;
elseif iscell(filenames)
   numf = size(filenames,2);
   dt = 2;
   filenames = sort(filenames); % this needs to happen; not sure why Matlab doesn't read in files in the same alphabetical order as selected in the UI

end

% for pre-processed Gait data, how many subjects do we have?

if db == 11
    file = [pathname filenames];
    filedata = load(file);

    heel_all = filedata.Heel;
    numh = size(heel_all,1);
   % swing_pdiff_all = filedata.swing_pdiff; % the signed value: positive means L > R; this is a cell array
    
    numf = size(heel_all,1);
    
    toe_all = filedata.Toe;
    numt = size(toe_all,1);
    
    gender_all = filedata.Gender; % 1 = male, 2 = female
    age_all = filedata.Age;
    
    % check heel and toe
    if numh == numt
        % OK
    else
        fprintf('\n Warning: number of cases of Heel and Toe data do not match!\n\n')
        return
    end
       
end

cd(pathname)

fprintf(['\n Read in ' num2str(numf) ' files.\n'])

% do we want to auto-exclude outliers identified by Hausdorff or RJE measure?
ex_method = input('\n Outlier exclusion method: \n   [0] Load existing excluded indices \n   [1] Hausdorff auto (3SD > median) \n   [2] Hausdorff auto + manual \n   [3] Ellis auto (successive ratio) \n   [4] Ellis auto + manual <ENTER> \n   [5] None \n   --> ');

if isempty(ex_method)
    ex_method = 4;
end

if ex_method == 0
    % assumes we have this as an available matrix
    [exname, expathname] = uigetfile( {filt,'.mat file'},'Select the matrix of excluded timestamps','MultiSelect', 'off');
    ex_file  = [expathname exname];
    ex_data = load(ex_file);
    ex_cap  = ex_data.ex_cap; % so there is no confusion about indexing!
    ex_inds = ex_data.ex_inds;
    
    % make sure it is the right size
    if size(ex_inds,1) == numf
        % OK
    else
        fprintf('\n *** Error: the number of subjects in the "excluded indices" file does not match the number of data files.\n     Respecify. \n\n')
        return
    end
else
    ex_inds = cell(numf,1); % to save the start and stop indices
    ex_cap = input('\n Number of timestamps to exclude at the beginning and end of file (<ENTER> for 4): ');
end

if isempty(ex_cap)
    ex_cap = 4;
end

% We may just want to look at individual files (even after we've specified
% the excluded time stamp matrix)

select_f = input('\n Select individual file(s) for examination? \n   [1] Yes \n   [2] No <ENTER> \n   --> ');

if isempty(select_f)
    select_f = 2;
end

if select_f == 1
    cases = input('\n Enter the set of cases to examine, in [ ]: ');
    numf = numel(cases);
else
    % use all cases
    cases = 1:numf;
end

  
%% Other parameters

% how large is the deletion window?
if ex_method ~= 0
    del_itis = input(' How many ITIs around each point to delete? (0 to 4) <ENTER uses 0>: ');
elseif ex_method == 0
    del_itis = 0; % we already marked the cases
end

if isempty(del_itis)
    del_itis = 0;
end

% downsample?
ds_choice = input(' Downsample the timeseries? \n   [1] Yes \n   [2] No <ENTER> \n   --> ');

if isempty(ds_choice)
    ds_choice = 2;
end

if ds_choice == 1
    HZ = input(' Enter the target sampling rate, in Hz: ');
    fprintf('\n') % leave this here
else
    HZ = Inf;
end

% plot?
if ex_method == 2 || ex_method == 4 % have to plot for safety's sake
    make_plots = 1;
else
    make_plots = input(' Plot graphs and pause between files? \n   [1] Yes \n   [2] No <ENTER> \n   --> ');
end

if isempty(make_plots)
    make_plots = 2;
end

% how long to pause
if make_plots == 1
    pause_dur = input('\n Duration to pause between viewing files, in seconds (<ENTER> for 1.0): ');
    
    if isempty(pause_dur)
        pause_dur = 1.0;
    end
end

% do DFA?
do_dfa = input(' Perform DFA analysis? \n   [1] Yes <ENTER> \n   [2] No \n   --> ');

if isempty(do_dfa)
    do_dfa = 1;
end

% only look at first N steps
num_events = input(' Exact number of steps to use (<ENTER> uses all): ');

if isempty(num_events)
    num_events = Inf;
end

if num_events == Inf
    % exclude subjects that don't have min number of steps
    min_events = input(' Minimum number of steps to retain a subject (<ENTER> retains all subjects): ');

    if isempty(min_events)
        min_step_thr = 0; % any positive number is fine
    end

else
    min_step_thr = 0;
end

%new_name = input(' Enter a character string descriptive name, or <ENTER> for none: ','s');

%if isempty(new_name);
%    use_name = [];
%else
    use_name = ['su' num2str(numf) '_hz' num2str(HZ) '_ev' num2str(num_events)];
%end

PRC = 50:10:90; % default is 50:10:90
    
% variables defined

    file_plus_times = cell(numf,3); % name, R, L
    file_name   = cell(numf,1);
    samp_time_min = nan(numf,1); % the *observed* time between sampled points
    samp_time_med = nan(numf,1);
    samp_time_max = nan(numf,1);
    
    % previously used stats
    walk_time   = nan(numf,1);
    walk_percent = nan(numf,1);
    step_count  = nan(numf,1);       
    cadence     = nan(numf,1); 
    pciL        = nan(numf,1); % Hausdorff
    pciR        = nan(numf,1); % Hausdorff
    pciP        = nan(numf,1); % the "correct" one, based on swing time calcs (per Plotnik)
    pciMax        = nan(numf,1); % the higher of the two values (not merged data)
    pciMer        = nan(numf,1); % merged L and R data

    cvsd       = nan(numf,1); 
    cvrms      = nan(numf,1);

    dfa_alpha = nan(numf,1);  
    
    % new, by RJE
    loc_prc      = nan(numf,numel(PRC));
    loc_transf_mn = nan(numf,1);
    lm_mean_loc  = nan(numf,numel(PRC)); 
    um_mean_loc  = nan(numf,numel(PRC)); 
    dm_mean_loc  = nan(numf,numel(PRC)); 
    
    suc_prc      = nan(numf,numel(PRC));
    suc_transf_mn = nan(numf,1);
    lm_mean_suc  = nan(numf,numel(PRC)); 
    um_mean_suc  = nan(numf,numel(PRC)); 
    dm_mean_suc  = nan(numf,numel(PRC)); 
    
    % swing time
    swingL_mn = nan(numf,1);
    swingL_cv = nan(numf,1);
    swingR_mn = nan(numf,1);
    swingR_cv = nan(numf,1);
    swingMer_mn = nan(numf,1);
    swingMer_cv = nan(numf,1);
    swingMax_mn = nan(numf,1);
    swingMax_cv = nan(numf,1);
    
    % counting up outliers  
    num_IEI              = nan(numf,1);
    inc_prc_H            = nan(numf,1); % Hausdorrf method   
    inc_prc_E            = nan(numf,1); % RJE method  

    % to save the group code; default = 0, which will ignore that subject
    
    group_code = zeros(numf,1);

%% loop through all files

fprintf(' Working ...\n')

% do this for physionet extractions

progressbar % Init single bar

for f = 1:numf
    casef = cases(f);
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if db <= 10 % for all methods

            progressbar(f/numf) % Update progress bar

            % open the file
           if dt == 1
               fname = filenames(casef,:);
           elseif dt == 2
               fname = char(filenames(casef));
           end

           % now we want to code our Experimental (= 2) and Control (= 1) subjects
           if sum(findstr(fname,grpC)) > 0
               group_code(f) = 1; % 1 = HC
           elseif sum(findstr(fname,grpE)) > 0
               group_code(f) = 2; % 2 = PD
           end

           file = [pathname fname]; 

           filedata = load(file);
           
           % *************************
           % now specific methods
           if db == 1
                % sampling rate is 100 Hz
                time = filedata(:,1);
                
                lh = filedata(:,2); % left heel
                rh = filedata(:,10); % right heel
                
                % sampling intervals
                dtime = diff(time);
                plot_dtime = 1;
                if plot_dtime == 1
                    figure(400)
                    plot(sort(dtime))
                    pause(1)
                    clf
                end
                
                samp_time_min(f) = min(diff(time));
                samp_time_med(f) = median(diff(time));
                samp_time_max(f) = max(diff(time));

                seg_full = [time lh rh];

                % rje: needs to be integrated with Hang's code as of 2013-08-22
                ex = [];
                [Xmat_old L R] = strike_extract(seg_full,1,fname,make_plots,ex); % L and R are used by td_outlier_detect
           end
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif db == 11
        casef = cases(f);
        % loop through all cases

        progressbar(f/numf) % Update progress bar

        fname = heel_all{casef,1}; % this is the file name

        % now we want to code our Experimental (= 2) and Control (= 1) subjects
           if sum(findstr(fname,grpC)) > 0
               group_code(f) = 1;
           elseif sum(findstr(fname,grpE)) > 0
               group_code(f) = 2;
           end

        % get the L and R (heel) steps

        L = heel_all{casef,2};
        R = heel_all{casef,3};
        
        % get L and R toes
        Lt = toe_all{casef,2};
        Rt = toe_all{casef,3};
        
        toe_off.L = Lt;
        toe_off.R = Rt;
        
        % update 2013-08-11: we work with the stride time series first to
        % do outlier clean up, and then 
        
       
        
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif db < 20
       % yet to be written
       
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    elseif db < 30
       % yet to be written 
       
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TEXT FILE METHODS
    
    elseif db >= 91
        % open the file
        progressbar(f/numf) % Update progress bar
        if dt == 1
           fname = filenames(casef,:);
        elseif dt == 2
           fname = char(filenames(casef));
        end
        
        filedata = load(fname);

        % are we in seconds or ms?
        if median(filedata) > 100
            % milliseconds
            filedata = filedata / 1000;
        else
            % seconds; OK
        end
        
        % if we have IEIs, we need to get back to timestamps for outlier
        % detection to work properly
        
        if db == 91
            % already timestamps
            L = filedata;
            R = [];
        elseif db == 92
            % need to get timestamps
            L = cumsum(filedata);
            R = [];
        end
        
    end % END OF ALL PROCESSING METHODS

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % OUTLIER DETECTION

    if ex_method == 0
        ex_inds_import = ex_inds{casef};
    else
        ex_inds_import = []; % 
    end

    long_run_method = 'c'; % concatenate runs for Physionet, but not for BEATS
    run_dur_thr = 10; % in seconds
    gap_dur_thr = Inf; % we don't care about distance between gaps
    loc_meth = 'kde_botev'; 
    outlier_meth = 'pc'; % percent change
    outlier_thr = 50; % 
    do_slope = 0; % don't need it
    % rje: we can probably make this number lower after looking at real data!
    [Xmat stats] = td_outlier_detect(L, R, ex_cap, ex_method, ex_inds_import, del_itis, long_run_method, loc_meth, outlier_meth, outlier_thr, run_dur_thr, gap_dur_thr, do_slope, make_plots,casef,fname); % Xmat will get fed DIRECTLY into td_calcs
    % note: long_run_ts will run inside of the above script

    % keep the excluded indices
    ex_inds{f} = stats.ex_inds; % note: in fact, we can load in the full matrix of ex_inds_import, 
                                % but then *replace* certain cases automatically when we write it *back* in this step!

    % how many outliers per the Hausdorff vs RJE method (not including
    % other manual exclusions, left-right alternations, long runs, etc)
    inc_prc_H(f) = 100* sum(Xmat(:,3))/numel(Xmat(:,3)); 
    inc_prc_E(f) = 100* sum(Xmat(:,4))/numel(Xmat(:,4)); % sum the 1s to get number of events in all the good runs

    % OK, now we still have to just pull out the RUNS; loop through each run separately
    % code has been modified to do this; will identify the runs in col4
    % and extract the data, and do that for all runs

    % 2013-09-05
    % create the *downsampled* version of the real data according to user inputs
    
    Xorig = Xmat(:,1); % timestamps
    
    if ds_choice == 2
        % don't downsample
        Xds = Xorig;
    elseif ds_choice == 1
        sr = (1 / HZ); % samples every xx *seconds*
        Xds = ceil(Xorig / sr) * sr; % the new time stamp is rounded *up* to the next sampling point
    end
    
%     % compare original vs downsampled
%     figure(41)
%     plot(diff(Xorig),'b')
%     hold on
%     plot(diff(Xds),'r')
%     hold off
    
    % now put it back in the matrix; we only test the downsampled version (if at all)   
    Xmat(:,1) = Xds;

    % now we do calculations after we have removed outliers, checked for alternating events (if relevant)
    
    %%%%%%%%%%%%%%%%%%%%%%%
    % PCI calculations now

    %cur_file = f
    plot_phi = 0; % 0 = don't plot; 1 = plot
    cur_file = f;
    toe_off = []; % will ignore this for now
    pci_res = pci_calc_concat(Xmat,plot_phi,num_events, toe_off); % will calculate both pciL and pciR; we don't need to select which one is correct within this step; do it later

    %%%%%%%%%%%%%%%%%%%%%%%
    % remaining td calculations

    location_meth = 'tm'; % trimean
    transform_meth = 'pc'; 

    plot_td = make_plots;
    use_interp = 0; %2013.10.28 - don't need to do this
    fignum = 150;
    td_res = td_calcs_concat(Xmat, PRC, location_meth, transform_meth, do_dfa, num_events, plot_td, use_interp,fignum); % will do calcs on all stable segments within the full time series
    
    if make_plots == 1
        pause(pause_dur) % leave this here; it is not an input parameter to td_calcs_concat, for simplicity
    end
        
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STORE VARIABLES

    file_name{f}         = fname;
    file_plus_times{f,1} = fname;
    file_plus_times{f,2} = R;
    file_plus_times{f,3} = L;
    
    pciL(f)             = pci_res.pciL;
    pciR(f)             = pci_res.pciR;
    
    swingL_mn(f)        = pci_res.swingL_mn;
    swingL_cv(f)        = pci_res.swingL_cv;
    swingR_mn(f)        = pci_res.swingR_mn;
    swingR_cv(f)        = pci_res.swingR_cv;
    swingMer_mn(f)         = pci_res.swingMer_mn;
    swingMer_cv(f)         = pci_res.swingMer_cv;
    swingMax_mn(f)         = pci_res.swingMax_mn;
    swingMax_cv(f)         = pci_res.swingMax_cv;
    
    
    %% now we determine whether pciL or pciR is the "correct" one per Plotnik method
    
    % 2013.10.28: The method below only works if we import the Heel + Toe +
    % Swing data right now
    
    if db == 11
%         if swing_pdiff_all{f} >= 0
%             pciF(f) = pciL(f);
%         else
%             pciF(f) = pciR(f);
%         end

        if swingL_mn(f) >= swingR_mn(f)
            pciP(f) = pciL(f);
        else
            pciP(f) = pciR(f);
        end

        
        % now just take the higher of the two values (i.e., in case we want
        % to use PCI but don't have swing time data)
        pciMax(f) = max(pciL(f),pciR(f));
    else
        pciP(f) = NaN; % irrelevant
        pciMax(f) = NaN;
    end
    
    pciMer(f)             = pci_res.pciMer;

    num_IEI(f)          = td_res.num_total; % how many total events
    step_count(f)       = td_res.num_incl; % tallys the number of retained (good) events
    walk_time(f)        = td_res.duration; % just on good data
    walk_percent(f)     = td_res.percent_valid; % percent of total duration of walk
    cadence(f)          = 60 / td_res.trimean; % for Physionet, location can be trimean
    
    cvsd(f)             = td_res.cvsd;
    cvrms(f)            = td_res.cvrms;
    dfa_alpha(f)        = td_res.dfa_alpha;

    loc_prc(f,:)        = td_res.loc_prc;
    loc_transf_mn(f)    = td_res.loc_transf_mn;
    lm_mean_loc(f,:)    = td_res.lm_mean_loc;
    um_mean_loc(f,:)    = td_res.um_mean_loc;
    dm_mean_loc(f,:)    = td_res.dm_mean_loc;
    
    suc_prc(f,:)        = td_res.suc_prc;
    suc_transf_mn(f)    = td_res.suc_transf_mn;
    lm_mean_suc(f,:)    = td_res.lm_mean_suc;
    um_mean_suc(f,:)    = td_res.um_mean_suc;
    dm_mean_suc(f,:)    = td_res.dm_mean_suc;
    
    
    
    
end % end of file loop
progressbar(1) % make sure it is closed

%% output variables

% compare the H and E methods for retained events
inc_prc = [inc_prc_H inc_prc_E (inc_prc_H-inc_prc_E)];

% finally we need to output all this as a single matrix for mc_2samp to run
% 2013-09-17: only use the "official" pci measure for simplicity
outcomes = [cadence cvsd cvrms swingMer_mn swingMer_cv swingMax_mn swingMax_cv pciP pciMax dfa_alpha loc_transf_mn suc_transf_mn];
other = [step_count walk_time walk_percent]; % separate out non-diagnostic features

% a few final exclusions
group_code(step_count < min_step_thr) = NaN; % if step_count = 0, then this will also return NaN


% just put three things into the output
output.down_samp_rate   = HZ;
output.num_events       = num_events;
output.file_name        = file_name;

output.file_plus_times  = file_plus_times;
output.samp_time_min    = samp_time_min;
output.samp_time_med    = samp_time_med;
output.samp_time_max    = samp_time_max;
output.group_code       = group_code;
output.gender           = gender_all;
output.age              = age_all;
output.outcomes         = outcomes;
output.other            = other;
output.num_IEI          = num_IEI;
output.inc_prc          = inc_prc;

ex_output.file_name = file_name;
ex_output.ex_cap    = ex_cap;
ex_output.ex_inds   = ex_inds;

% save date and time
save_date = datestr(now,29); % the date as YYYY-MM-DD
save_time = datestr(now,13); % the time
save_time = [save_time(1:2) '.' save_time(4:5) '.' save_time(7:8)]; % so it will be a saveable name
save_dt = [save_date '_' save_time];


% if isempty(new_name)
%     save_inst1 = ['save physionet_results_' save_dt '.mat output'];
%     save_name1 = [     'physionet_results_' save_dt '.mat'];
%     
%     save_inst2 = ['save physionet_ex_inds_' save_dt '.mat ex_output'];
%     %save_name2 = [     'physionet_ex_inds_' save_dt '.mat'];
%     
% else
    save_inst1 = ['save ' use_name 'physionet_results_' save_dt '.mat output'];
    save_name1 = [        use_name 'physionet_results_' save_dt '.mat'];

    save_inst2 = ['save ' use_name 'physionet_ex_inds_' save_dt '.mat ex_output'];
    %save_name2 = [        use_name 'physionet_ex_inds_' save_dt '.mat'];    
% end

eval(save_inst1)  % that simple!
eval(save_inst2)


fprintf(['\n Finished. The *.mat output file was written to the working directory (pwd): \n <' save_name1 '> \n\n']);
