function output = cent_tend_brute(X,thr,min_prc,plot_it)

% brute-force measure of central tendency: find the value from which the
% most number of values have a percentage deviation of <= thr
%
% this only will work for POSITIVE vectors of numbers (e.g., inter-event intervals)
%
% "min_prc" is the minimum percentage of the series which constitutes a
% valid location; default = 0 (i.e., no minimum)
%
% rje | 2013.09.10

% timer start
tic;

X = X(:);
ind = 1:numel(X);

if nargin < 3
    min_prc = 0;
end

if nargin < 4
   plot_it = 0;
end

%% sort X - important for later

X = sort(X);

%% successive percent change (matrix version
X_o = X(1:end-1);
X_s = X(2:end);

suc_prc = 100 * (X_s - X_o) ./ X_o; % element wise division so we need ./

suc_prc = [NaN; suc_prc]; % to get indexing correct

%% find the longest run of gradual change

bin_suc_prc = suc_prc < thr; % 1 means good

outlier_meth = 'pc';
outlier_thr = Inf; % ignores it
run_dur_thr = -Inf; % no minimum length requirement
gap_dur_thr = 0; % no gaps allowed

long_run_out = long_run_ts(bin_suc_prc, outlier_meth, outlier_thr, run_dur_thr, gap_dur_thr);

% we only want the longest run
long_run_ind = long_run_out.long_run_ind;

%% safety measure
% if the number of elements in the longest run isnt > P
% percent of the data, then we use the median
long_run_prc = 100 * sum(long_run_ind) / numel(X);

prc_thr = min_prc;

if long_run_prc <= prc_thr % also works if prc_thr == 0
    location = median(X);

else
    % we are OK
    X_seg = X(long_run_ind == 1); % actual X-values in the longest run
    location = median(X_seg); % this way, the location is in the middle of the segment
end



%% plots
if plot_it == 1
   figure(500)
   subplot(2,1,1)
   plot(sort(X))
   hold on
   
   if long_run_prc > prc_thr
       % show the longest run
       plot(ind(long_run_ind == 1),X_seg,'r.')

       % show the median value
       plot(ind(X == location),X(X == location),'g.')
   end
   hold off
   
   subplot(2,1,2)
   plot(suc_prc)
   
   %figure(501)
   %plot(suc_prc,gradient(X),'.')
end

%% output

output.median = median(X);
output.mode = mode(X);
output.brute_location = location;
output.long_run_prc = long_run_prc;
output.unique_percent = 100 * numel(unique(X)) / numel(X); % how much repeated values do we have?

el_time = toc;

output.elapsed_time = el_time;


