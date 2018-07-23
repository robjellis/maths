function output = long_run_ts(Xmat, run_merge_meth, run_merge_thr, string_type, run_dur_thr, gap_dur_thr, do_slope)

%
% [indices] = long_run_ts(Xmat, run_dur_thr, gap_length)
% find the longest run of 1s in a binarized series (1s and 0s)
%
% note: long_run_ts.m is optimized for the Physionet project, not the BEATS
% project (which can work fine with IBIs)
%
% Xmat           = the output of td_outlier_detect.m    
% run_merge_thr  = depends on run_merge_meth [Inf will ignore this parameter]
% run_dur_thr    = minimum duration (in seconds) that a run should last [-Inf to ignore this param]
% gap_dur_thr    = maximum duration (in seconds) that a gap between runs may be [Inf to ignore this param]
%
% "output" will return a matrix of N rows x M columns; each column is a *distinct* run
% - unlike the BEATS project, we don't care how long gaps (e.g., turns) are
% 
% code by rje | version = 2014.03.27
% borrowing from
% http://stackoverflow.com/questions/3274043/finding-islands-of-zeros-in-a-sequence
%

%% start

% ***************************
% 2013.10.12: slope takes a long time to do on individual points, so there
% is now a method to do it on time-based sections

if nargin < 7
    do_slope = 0;
end

slope_meth = 'sec'; 
plot_slope = 0;

% ***************************

if run_dur_thr == 0 % special case; if we are using the entire excerpt, then run_dur_thr is set to 0
    param = 10;
else
    param = min(10,run_dur_thr); % we can't have a window longer than run_dur_thr
end


if do_slope == 1
    % check for robust regression (requires stats toolbox)
    checkrob = which('robustfit'); % make sure we have it

    if sum(size(checkrob))
        fit_method = 'r';
    else
        fit_method = 'p'; % polyfit (standard in Matlab)
    end
    
    % new 2014.07.02 - just use regular linear regression for simplicity
    use_regular = 1;
    
    if use_regular == 1
        fit_method = 'p';
    else
        % just leave fit_method as it was
    end

    % turn these warnings off; note: no regression will be calculated (slope_z will stay at NaN) 
    warning('off','stats:statrobustfit:IterationLimit')
    warning('off','stats:robustfit:RankDeficient')

    % number of points in the slope window
    pts = 10;
else
    do_slope = 0; % easy way to just cancel this operation
end

% we don't need this for physionet data but can leave it here anyway; just
% use Inf to ignore this value

if isempty(gap_dur_thr) % for Physionet
    gap_dur_thr = Inf; 
end



%% Next, find the starting indices, ending indices, and duration of each string of ones using the functions DIFF and FIND
% 1 = good, 0 = bad; these are for TIMESTAMP indices

% the assumption is that the *final* column of Xmat will have the 1s and 0s
% for good and bad indices

% if we just input a single vector of 1s and 0s, we need to do a little
% work here

if numel(unique(Xmat)) == 2
    % just 0s and 1s
    bin_data = Xmat;
    
    num_chan = 1;
    
    Xts = 1:numel(Xmat);
    Xts = Xts(:); % column
else
    
    bin_data = Xmat(:,end);

    chan_data = Xmat(:,2); % L = 1, R = 2

    % do we have two channels?
    num_chan = numel(unique(chan_data));
    
    Xts = Xmat(:,1); % timestamps

end

if num_chan == 1
    % then we don't care about alternations
    alt_chan = ones(size(Xmat(:,1))); % so we keep all indices
else  
    alt_chan = zeros(size(chan_data)); % to store valid alternations
end

dsig = diff(bin_data);
dsig = [0; dsig]; % so the indexes line up!


% be careful: we need to handle the case of [1 1 1 1 1 1 1 1 0 1 1 1 1 1 1 ... ]
if bin_data(1) == 1;
    startIndex = [1; find(dsig == 1)]; % so we keep the initial case!
else
    startIndex = find(dsig == 1); % find ALL changes from 0 to 1
end

endIndex = find(dsig == -1); % final ALL changes from 1 to 0
    
    if isempty(endIndex)
        endIndex(1) = numel(bin_data) + 1; 
    end
    
    if bin_data(end) == 1
        endIndex(end+1) = numel(bin_data) + 1; 
    end
    
    % make sure we don't have duplicate values
    endIndex = unique(endIndex);
    
    % make sure same vector orientation
    startIndex = startIndex(:);
    endIndex = endIndex(:);

    
endIndex = endIndex - 1; % the "-1" makes this work correctly since we want the *last* 1 *before* a 0

nruns = numel(startIndex);

run_num = (endIndex-startIndex) + 1; % will count number of elements in each successive run 

run_med_all = nan(1,nruns);
run_len_all = nan(1,nruns); % number of elements per run
run_dur_all = nan(1,nruns); % in *seconds*

% get the median value from within each run, using orig_data

for i = 1:nruns
   si = startIndex(i);
   ei = endIndex(i);
   seg =  Xts(si:ei);
   run_med_all(i) = median(diff(seg)); % median of inter-event interval
   
   % defined in seconds or events, depending on run_type
   if strcmp(string_type,'s')   
       run_dur_all(i) = seg(end) - seg(1); % in seconds
   elseif strcmp(string_type,'e') % duration defined by number of events
       run_dur_all(i) = numel(seg);
   end
   
   
   if num_chan == 2
       
       % what if si = ei?
       if si == numel(chan_data)
           % just leave this at zero
       else
           % also see if we have alternations within the runs
               if chan_data(si+1) ~= chan_data(si)
                   alt_chan(si:si+1) = 1; % both are 1
               end

           % now check the rest
           for a = si+2:ei    
               if chan_data(a) ~= chan_data(a-1)
                   % good
                   alt_chan(a) = 1; % keep this index
               else
                   % stay at zero
               end     
           end
       end % if si == ei
   end % num_chan == 2
end

%all_runs = []; % variable needs to be defined regardless of whether it is filled

%run_len

%% Then, find the strings of zeroes with a run_num greater than or equal to some value

% rje change: we only want to take long runs

stringIndex = run_dur_all >= run_dur_thr;  % this returns either 0 or 1
                                      % there may be mutiple excerpts that are long
                                      % this value effectively sets how widely
                                      % spaced beat "jitter" should be
nruns = sum(stringIndex);

% initialze the output vars - use ZEROS not NaNs
all_runs_ind = zeros(size(Xmat,1),1); % RJE: leave this as "zeros", not NaNs
long_run_ind = zeros(size(bin_data));
long_run_dur = NaN(size(bin_data));
gap_dur_save = 0; % needs to be zero to make life easier (calculating stable percentage, etc.)
slope_pc_max_abs = NaN;


% are we in luck?

% initialize variables
keep_runs = []; % will tell us which runs (e.g., [2 3 4]) we end up keeping
slope_pc_max_abs = [];
x_start = [];
x_end = [];
y_start = [];
y_end = [];
    
if sum(stringIndex) == 0
    % output indicies stay at zero
    long_run_dur = 0;
elseif sum(stringIndex) >= 1
    % how many do we have?
    
    
    % just keep the indices that satisfy run_dur_thr

    startIndex  = startIndex(stringIndex);
    endIndex    = endIndex(stringIndex);
    run_med_keep = run_med_all(stringIndex);
    run_dur_keep = run_dur_all(stringIndex);
    
    
    % now find the gaps between the *retained* runs
    
    if nruns > 1
        for i = 1:nruns-1
           gap_dur_keep(i) = Xts(startIndex(i+1)) - Xts(endIndex(i));
        end
    else
        gap_dur_keep(1) = NaN; % only one Run, so no Gaps
    end
    
    startIndex_orig = startIndex;
    endIndex_orig = endIndex;
    

    
    if nruns == 1
        % to save time we can do this in one step
        all_runs_ind(startIndex:endIndex) = 1;
        long_run_ind = all_runs_ind;
        long_run_dur = Xts(endIndex) - Xts(startIndex);
        keep_runs = 1;
        
        if do_slope == 1
            run = Xts(startIndex:endIndex);
            drun = diff(run);

            % slope
            
            slope_out = robustfit_slope(run(2:end),drun,slope_meth,param,fit_method,plot_slope);
            slope_pc_max_abs = slope_out.slope_pc_max_abs;
            x_start = [x_start; slope_out.x_start];
            x_end = [x_end; slope_out.x_end];
            y_start = [y_start; slope_out.y_start];
            y_end = [y_end; slope_out.y_end];
        else
            slope_pc_max_abs = NaN;
        end
        
    elseif nruns > 1
        % save the locations of these runs BEFORE we worry about gaps
        for nn = 1:nruns
            all_runs_ind(startIndex(nn):endIndex(nn)) = 1; % will be a column  
        end
        
        % Now, look for consecutive long runs that may be separated by a max of gap_length events

        %all_runs = xdiff(startIndex(1):endIndex(1)); % the first run
        if do_slope == 1
            run = Xts(startIndex_orig(1):endIndex_orig(1)); % the first run
            drun = diff(run);
            % slope
            slope_out = robustfit_slope(run(2:end),drun,slope_meth,param,fit_method,plot_slope);
            
            slope_pc_max_abs(1) = slope_out.slope_pc_max_abs; % will increment this
            x_start = [x_start; slope_out.x_start];
            x_end = [x_end; slope_out.x_end];
            y_start = [y_start; slope_out.y_start];
            y_end = [y_end; slope_out.y_end];            
        else
            slope_pc_max_abs = NaN;
        end
        
        for p = 2:nruns
            
            ei = endIndex(p-1);
            si = startIndex(p); 
            
            % switch between seconds and number of events
            if strcmp(string_type,'s')
                gap_dur = Xts(si) - Xts(ei); % time between the previous run's last good beat and the next run's first good beat
            elseif strcmp(string_type,'e')
                gap_dur = numel(ei+1:si-1);  % number of elements; % the plus 1 and minus 1 are key since we want the time *between* Run(i) and Run(i+1)           
            end
            % 2013-08-23: defined as a ratio, just like do for beats
            if strcmp(run_merge_meth,'q')
                run_quot = max(run_med(p),run_med(p-1)) / min(run_med(p),run_med(p-1)); % divde larger by smaller to get a quotient >= 1.0
            
                if gap_dur <= gap_dur_thr && run_quot <= run_merge_thr % make sure successive runs have (1) tolerable spacing and (2) tolerable difference in median
                    % we now remove this gap
                    endIndex(p-1) = NaN;
                    startIndex(p) = NaN;
                    gap_dur_save(p) = gap_dur; % create a running list of the time of *included* gaps 
                    % since we know the gap is OK, we can now add the NEXT run
                    
                    %next_run = xdiff(startIndex_orig(p):endIndex_orig(p));
                    %all_runs = [all_runs; next_run]; % a running column; now we are just keeping the good data
                    
                    if do_slope == 1
                        run = Xts(startIndex_orig(p):endIndex_orig(p)); % the next run
                        drun = diff(run);
                        
                        % slope
                        slope_out = robustfit_slope(run(2:end),drun,slope_meth,param,fit_method,plot_slope);
                        slope_pc_max_abs(p) = slope_out.slope_pc_max_abs; % will increment this 
                        x_start = [x_start; slope_out.x_start];
                        x_end = [x_end; slope_out.x_end];
                        y_start = [y_start; slope_out.y_start];
                        y_end = [y_end; slope_out.y_end];                        
                    else
                        slope_pc_max_abs(p) = NaN;
                    end
                    
                else
                    % do nothing
                end
            elseif strcmp(run_merge_meth,'pc') % percent change
                run_pc = 100 * (run_med(p) - run_med(p-1)) / run_med(p-1); % percent difference
            
                if gap_dur <= gap_dur_thr && run_pc <= run_merge_thr % make sure successive runs have (1) tolerable spacing and (2) tolerable difference in median
                    % we now remove this gap
                    endIndex(p-1) = NaN;
                    startIndex(p) = NaN;
                    gap_dur_save(p) = gap_dur; % create a running list of the time of *included* gaps   
                    % since we know the gap is OK, we can now add the NEXT run
                    
                    %next_run = xdiff(startIndex_orig(p):endIndex_orig(p));
                    %all_runs = [all_runs; next_run]; % a running column; now we are just keeping the good data
                    
                    if do_slope == 1
                        run = Xts(startIndex_orig(p):endIndex_orig(p)); % the next run
                        drun = diff(run);

                        % slope
                        slope_out = robustfit_slope(run(2:end),drun,slope_meth,param,fit_method,plot_slope);
                        slope_pc_max_abs(p) = slope_out.slope_pc_max_abs; % will increment this 
                        x_start = [x_start; slope_out.x_start];
                        x_end = [x_end; slope_out.x_end];
                        y_start = [y_start; slope_out.y_start];
                        y_end = [y_end; slope_out.y_end];                        
                    else
                        slope_pc_max_abs(p) = NaN;
                    end
                    
                else
                    % do nothing
                end 
                
            elseif strcmp(run_merge_meth,'pd')
                run_pd = 100* (run_med(p) - run_med(p-1)) / (0.5 * (run_med(p) + run_med(p-1))); % percent difference
            
                if gap_dur <= gap_dur_thr && run_pd <= run_merge_thr % make sure successive runs have (1) tolerable spacing and (2) tolerable difference in median
                    % we now remove this gap
                    endIndex(p-1) = NaN;
                    startIndex(p) = NaN;
                    gap_dur_save(p) = gap_dur; % create a running list of the time of *included* gaps   
                    % since we know the gap is OK, we can now add the NEXT run
                    
                    %next_run = xdiff(startIndex_orig(p):endIndex_orig(p));
                    %all_runs = [all_runs; next_run]; % a running column; now we are just keeping the good data
                else
                    % do nothing
                end 
            elseif strcmp(run_merge_meth,'free')
                % we don't care about differences between runs; we only
                % care about gap_dur
                
                if gap_dur <= gap_dur_thr 
                    % we now remove this gap
                    endIndex(p-1) = NaN;
                    startIndex(p) = NaN;
                    gap_dur_save(p) = gap_dur; % create a running list of the time of *included* gaps   
                    
                    if do_slope == 1
                        run = Xts(startIndex_orig(p):endIndex_orig(p)); % the next run
                        drun = diff(run);
                        
                        % slope
                        slope_out = robustfit_slope(run(2:end),drun,slope_meth,param,fit_method,plot_slope);
                        slope_pc_max_abs(p) = slope_out.slope_pc_max_abs; % will increment this  
                        x_start = [x_start; slope_out.x_start];
                        x_end = [x_end; slope_out.x_end];
                        y_start = [y_start; slope_out.y_start];
                        y_end = [y_end; slope_out.y_end];                        
                    else
                        slope_pc_max_abs(p) = NaN;
                    end                    
                else
                    % do nothing
                end 
            end

        end

        % get rid of NaN
        startIndex = startIndex(startIndex > 0);
        endIndex   = endIndex(endIndex > 0);

        
        % now we can get the final vector, with 1 = valid and 0 = invalid

            if numel(startIndex) == 0
                % stays at zero
            else

                % get the longest run (this is not necess. for Physionet at this point)

                diffInd = endIndex - startIndex;
                start_long = min(startIndex(diffInd == max(diffInd))); % we just want one excerpt
                end_long =   min(endIndex  (diffInd == max(diffInd)));  % we just want one excerpt       

                long_run_ind = zeros(size(bin_data));
                long_run_ind(start_long:end_long) = 1;
                
                long_run_dur = Xts(end_long) - Xts(start_long);
            end

    end
end

%% now incorporate the alternating index
% may need to work this in better to the above (2014.03.27)

all_runs_ind = min(alt_chan,all_runs_ind);
long_run_ind = min(alt_chan,long_run_ind);

%% outputs
output.all_runs_ind = all_runs_ind; % identifies the indices of all runs
output.long_run_ind = long_run_ind; % 1s and 0s at index locations for the *longest* run (after ignoring gaps)
output.long_run_dur = long_run_dur;
output.gap_dur_sum = sum(gap_dur_save);
output.slope_pc_max_abs = max(slope_pc_max_abs); % single value for the file

% group these for easy export
if do_slope == 1
    plot_points.x_start = x_start;
    plot_points.x_end = x_end;
    plot_points.y_start = y_start;
    plot_points.y_end = y_end;
else
    plot_points.x_start = NaN;
    plot_points.x_end = NaN;
    plot_points.y_start = NaN;
    plot_points.y_end = NaN;
end

output.plot_slope_points = plot_points;

