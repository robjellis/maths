function [output] = irace_accel

% analysis of triaxial accelerometer data from iPhone
%
% methods based on:
% - Zijlstra 2003, 2004
% - Yang et al 2012 (iGAIT)
%
% last update = 2013.04.24
% RJE

% select files
[filenames, pathname] = uigetfile( {'*.txt','TXT-files (*.txt)'},'Select file(s)','MultiSelect', 'on');

% how many files?
if ischar(filenames)
   numf = 1;
   dt = 1;
elseif iscell(filenames)
   numf = size(filenames,2);
   dt = 2;
end

% other input parameters

chan = input(' Plot channels? [1] Yes <ENTER>; [2] No>: ');
if isempty(chan)
    chan = 1;
end

posit = input(' Device position: [1] Back; [2] Front <ENTER>: ');
if isempty(posit)
   posit = 2;   
end

if posit == 1
    pos = -1; % multiply the channel by -1
elseif posit == 2
    pos = 1; % i.e., don't flip the channel
end


filter = input(' Filtering option: [1] Butterworth; [2] Sgolay; [3] moving average; [4] None <ENTER>: ');
if isempty(filter)
    filter = 4;
end


for f = 1:numf
   if dt == 1
       fname = filenames(f,:);
   elseif dt == 2
       fname = char(filenames(f));
   end
       
   file = [pathname fname]; 
    
   X = load(file); 
   
   % figure out if we are dealing with sync/cont or uncued
   
   if size(X,2) == 4
       % uncued
       shiftcol = 0;
   elseif size(X,2) == 5
       % sync/cont
       shiftcol = 1;
   end
       
   
   % get the sampling rate (in Hz) from the timestamp (ts) series
   ts = X(:,1+shiftcol);
   tsd = diff (ts);
   tsmax = max(ts);
   tscum = cumsum(ts); % get the cumulative sum of time
   
   % cut off the beginning and ends of the file (in seconds)
   cutsec = 1;
   tsmin = find(ts > cutsec, 1 );
   tsmax = find(ts < (tsmax - cutsec), 1, 'last' );
   
   % now get that range for X
   ts = ts(tsmin:tsmax);
   X = X(tsmin:tsmax,:);
      
   % plot the actual Hz rates
   %figure(10)
   %plot(1./tsd)
   
   tsmed = median(tsd); % the median time stamp - unbiased
   sr = 1/tsmed; % in Hz
   
   % get the antero-posterior acceleration
   % in iPhone, this is the 4th column of data
   
   lr = X(:,2+shiftcol); %lr = lr - median(lr); % left-right
   ud = X(:,3+shiftcol); %ud = ud - median(ud); % vertical
   ap = X(:,4+shiftcol);
   
   if chan == 1
      % rje confirmed these labels are correct
      figure(20)
      subplot(3,2,1)
      plot(ts,lr,'g')
      ylabel('L-R')
      
      subplot(3,2,3)
      plot(ts,ud,'r')
      ylabel('U-D')
      
      subplot(3,2,5)
      plot(ts,ap,'b')
      ylabel('A-P')
   end
       
   
   % multiply ap by -1 if iPod is positioned on the subject's back (spine)
   ap = ap * pos;
      
   norm = sqrt(ud.^2 + lr.^2 + ap.^2);
   
   % create the low-pass 4th-order Butterworth filter with cut off of 5 Hz
   % note that the cufoff frequency is NORMALIZED; see
   % http://www.mathworks.com/matlabcentral/newsreader/view_thread/263952
   
   % first get the range from -1 to 1
   
   ap = ap - min(ap); % min is 0
   ap = ap / max(ap); % max is 1
 
   
   if filter == 4
      % no filtering
      apfilt = ap;
      
   elseif filter == 1
       % set up the filter
       order = 4;
       cutoff = 5; % Hz
       [b,a]=butter(order,(2*pi*cutoff)/sr,'low');

       % apply the filter ... use filtfilt to get zero-order
       apfilt = filtfilt(b,a,ap);
       
   elseif filter == 2  
       apfilt = smooth(ap,9,'sgolay',3); 
   elseif filter == 3
       apfilt = smooth(ap,3); % simple moving average
   end
   
   subplot(3,2,[2 4 6])
   plot(ts,ap,'k',ts,apfilt,'m')
   xlabel('Time (s)')
   ylabel('A-P acceleration trace')
   
   % peak detection
   delta = 0.2;      % note that the full range of values is 0 to 1
   [maxtab, mintab] = peakdet(apfilt, delta,ts);
   hold on; plot(mintab(:,1), mintab(:,2), 'g*');
   plot(maxtab(:,1), maxtab(:,2), 'b*');
   
   
   stept = maxtab(:,1); % the timestamps, in seconds
   stepv = maxtab(:,2); % timestamp values (y-axis)
   
   % -----------------------------------------------
   % for stept 2 to N-1, we look on either side to see if the distance is small
   
   % initialization
   nstep = numel(stept);
   minT = .300; % in seconds; minimum time between sucessive steps
   stepind = 1:nstep; % just for the initial iteration
   g = 2;
   
   for n = 2:nstep-1
       look = g-1 : g; % use the real indices
       tt = stept(look);
       vv = stepv(look);
       
       if diff(tt) >= minT
           % we're fine; both are greater than minT
           g = g+1;
       else
           % exclude the lower voltage value
           
           ex = min(look(vv == min(vv))); % the index value to exclude; if there is a tie, we exclude the first one 
           keep = stepind(stepind ~= ex);
           
           stept = stept(keep);
           stepv = stepv(keep);
          
           % do not advance the counter
           
           stepind = 1:numel(stept); % reset this
       end
       
           if g > nstep
              break
           end
   end
   
   
   % now we plot the retained steps
   plot(stept, stepv, 'ro');
   hold off
   
   stept_diff = diff(stept);
   
   figure(21)
   plot(stept_diff,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0],'Marker','.')
   ylabel('Inter-strike Interval (s)')
   xlabel('Ordinal events')
   
   if numf > 1
        pause(1)
   end
   
   ISI{f} = stept_diff;
   
end % file loop

output.ISI = ISI;
