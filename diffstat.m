function output = diffstat(X, local_thresh, run_dur_thresh, gap_dur_thresh)

% an intelligent way of finding the central tendency in a time series that
% has a rhythmic process (i.e., an underlying base IOI)

if nargin < 2
    local_thresh = 5; % a liberal perceptual beat-by-beat benchmark
end

if nargin < 4
    run_dur_thresh = 10; % in seconds
    gap_dur_thresh = 2;  % in seconds
end

X = X(:);
nx = numel(X);

% updated on 19 July

Xadm = zeros(nx,1);
Xasd = zeros(nx,1);

Xadm(1) = 1; % assume that the first element could start a run
Xasd(1) = 1;

% note: this will not detect, for example, a gradual change over time; all
% successive differences may be sub-threshold. However, that is the point:
% we define LOCAL stability as separate from GLOBAL stability (PADM, PASD, etc.)

%Xmed = median(X);

% we DON'T want to use the median, because that may not actually represent
% the most frequent. Instead, we use the new central tendency statistic by
% RJE ("cent_tend.m")

ratio_thresh = 1/3; % default value

output = cent_tend(X,ratio_thresh);

loc = output.location;

for i = 2:nx;
   % percent absolute deviation from median
   Xadm(i)  = 100 * abs(X(i) - loc) / loc;
   
   % percent absolute successive difference 
   Xasd(i) = 100 * abs(X(i) - X(i-1)) ./ (0.5 * (X(i) + X(i-1)));
end

% now we do the binarizing option
   Xadm_bin = Xadm <= local_thresh;
   Xasd_bin = Xasd <= local_thresh;
   
   Xbin = (Xadm_bin + Xasd_bin) == 2;

out = long_run(X, Xbin, local_thresh, run_dur_thresh, gap_dur_thresh);

output.location = loc; 
output.indices = out.indices;
output.values = out.values;
output.duration = out.duration;
