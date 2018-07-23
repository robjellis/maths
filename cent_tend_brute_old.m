function output = cent_tend_brute(X,thr,minP,plot_it)

% brute-force measure of central tendency: find the value from which the
% most number of values have a percentage deviation of <= thr
%
% this only will work for POSITIVE vectors of numbers (e.g., inter-event intervals)
%
% rje | 2013.09.10

% timer start
tic;

X = X(:);
ind = 1:numel(X);

save_prc = zeros(numel(X),1);
suc_prc = nan(numel(X),1);

if nargin < 3
   minP = 0; % ignore this parameter
end

if nargin < 4
    plot_it = 0;
end

%% sort X - important for later

X = sort(X);

%% loop

for i = 1:numel(X)
   abs_prc = 100 * abs(X - X(i)) / X(i); 
   save_prc(i) = 100 * sum(abs_prc < thr) / numel(X);
   
end

for i = 2:numel(X)
    suc_prc(i) = 100 * (X(i) - X(i-1)) / X(i-1); % this is very threshold dependent!
end
%% make sure we have at least P percent of the data, otherwise use median

maxP = max(save_prc);

if maxP < minP
    % just use median
    location = median(X);
else
    max_inds = ind(save_prc == maxP);
    max_vals =   X(save_prc == maxP);
    
    location = max_vals(ceil(numel(max_inds)/2)); % this will give us an *actual* index value in case there is an even number of candidates; median won't work here
end



%% plots
if plot_it == 1
   figure(500)
   subplot(3,1,1)
   plot(sort(X))
   hold on
   plot(max_inds,max_vals,'r.')
   hold off
   
   subplot(3,1,2)
   plot(suc_prc)
   
   subplot(3,1,3)
   plot(save_prc)
   hold off

   
   %figure(501)
   %plot(suc_prc,gradient(X),'.')
end

el_time = toc;
%% output
output.location = location;
output.mode = mode(X);
output.median = median(X);
output.maxP = max(save_prc);
output.elapsed_time = el_time;

unique_prc = 100* numel(unique(X)) / numel(X)

