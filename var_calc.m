function [output] = var_calc(trim,ndev,win,sm_pr)

% function [output] = var_calc(trim,ndev,win)
%
% find the length of good event series for subsequent time-domain
% calculations
%
% trim = number of events to trim from the beg and end of each series
% ndev = outliers coded if > ndev median abs. deviations
% win = minimum length of a contiguous window free of artifacts
% sm_pr = value of the smoothing prior function; NaN will leave unsmoothed

% read in all the TXT files

[files path] = uigetfile({'*.txt','Text files (*.txt)'},'Select file(s)', 'MultiSelect', 'on');
%files = spm_select([1 Inf],'mat','Select file(s):',[],pwd,'.*');
numf = size(files,1);

size(files)

files
path
% file loop

tot_ev = nan(numf,1);
time1 = nan(numf,1);
time3 = nan(numf,1);
win3 = nan(numf,1);
num_bad = nan(numf,1);
num_win = nan(numf,1);
min_win = nan(numf,1);

mnITI = zeros(numf,1);
cvITI = zeros(numf,1);
cvrmsITI = zeros(numf,1);
dfaITI = zeros(numf,1);

% what sort of data format are we dealing with?
   file = char(files(1));
   file = [path file];
   %file = strtrim(files(1,:));  % strtrim removes white spaces
   data = load(file);
   tsize = size(data);
   % make sure correct size file!
   if tsize(2) == 1
       % then we are OK
       data = data(:);
       colch = 1;
   elseif tsize(2) > 1
       % then we have a matrix
       fprintf('\n File appears to be a matrix.\n')
       colch = input(' Enter the desired column to extract: ');
       
   end
   
for f = 1:numf
   % load the file
   file = char(files(f));
   file = [path file];
   %file = strtrim(files(f,:));  % strtrim removes white spaces
   data = load(file);
   fsize = size(data);
   
   % confirm correct size
   
   if tsize(2) ~= fsize(2) % check number of columns
      % bad
      fprintf('\n This file does not have the expected size. Respecify.\n\n');
      return
   else
       % OK
   
       % extract the desired column
       data = data(:,colch);

       % now we trim the file to ignore beginning and end
       datat = data(trim:end-trim+1);
       tot_ev(f) = numel(datat);

       % identify location of bad ITIs
       % note: mad is in the stats toolbox
       badev = (datat < (median(datat)-ndev*mad(datat))) + (datat > (median(datat)+ndev*mad(datat)));
       
       if sum(badev) == 0
           badev(1) = 1;
           badev(end) = 1; % just to help the code
       end

       % now, find the location of the first good window in the data  
       out = find_run(badev,win);

       win1 = out.start;
       win2 = out.stop;
       winL = out.length;
       win3(f) = winL;
       num_bad(f) = out.num_bad;
       num_win(f) = out.num_win;
       min_win(f) = out.min_win;

       if winL < win % if observed is less than target ...
          % stay at nan
       else
           % keep going

           dseg = datat(win1:win1+win-1); % just the ITIs in the window we want (from the trimmed data)

           % mean
           mnITI(f) = mean(dseg);

           % CV
           cvITI(f) = std(dseg) / mean(dseg) * 100;

           % CV RMS
           cvrmsITI(f) = rmssd(dseg) / mean(dseg) * 100;

           % DFA: [D,Alpha1] = DFA_main(DATA,n1,n2,steps) ... n1 to n2 is geometrically spaced
           [D alpha] = DFA_main(dseg,4,win,20); 
           dfaITI(f) = alpha;

           % now convert from events to actual time (seconds), including the trimmed part

           t1 = ceil(sum(data(1:(win1+trim-1))));
           t2 = floor(sum(data(1:(win1+win+trim-1)))); % note: this is to get the actual desired window size

           t3 = t2 - t1 + 1; % length of window in seconds
           % now convert to mins + secs
           t1m = floor(t1/60);
           t3m = floor(t3/60);

           t1s = rem(t1,60)/100; % remainder, in *decimal* form (i.e., .40 = 40 secs)
           t3s = rem(t3,60)/100;

           time1(f) = t1m + t1s;
           time3(f) = t3m + t3s;
       end
   
   end % fsize
end % f loop

output.tot_ev = tot_ev;
output.bad_ev = num_bad;
output.num_win = num_win;
output.min_win = min_win;
output.run_len = win3;

output.mnITI = mnITI;
output.cvITI = cvITI;
output.cvrmsITI = cvrmsITI;
output.dfaITI = dfaITI;

output.win_start = time1;
output.win_length = time3;

fprintf('\n Done. All window times are in mm.ss.\n\n');