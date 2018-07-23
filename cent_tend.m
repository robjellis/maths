function [output] = cent_tend(X,hist_ratio_thr)

% cent_tend(X,hist_ratio_thr)
% 
% A new central tendency (location) statistic (note: a "statistic" in a sample is an
% estimator of a "parameter" in a population). "X" is the inter-event
% interval series.
%
% This measure is a novel calculation; see
% http://en.wikipedia.org/wiki/Central_tendency for other measures
% See Matlab documention; both median and mode can be problematic in cases.
% Median:

%
% Mode: 
%    This function is most useful with discrete or coarsely rounded data.
%    The mode for a continuous probability distribution is defined as
%    the peak of its density function.  Applying the MODE function to a
%    sample from that distribution is unlikely to provide a good estimate
%    of the peak; it would be better to compute a histogram or density
%    estimate and calculate the peak of that estimate.  Also, the MODE
%    function is not suitable for finding peaks in distributions having
%    multiple modes.
%
% rje | 2013.03.28 

if nargin < 2
    hist_ratio_thr = 1/3; % that is, at least one third of of the time series values need to fit into this bin
end

X = X(:);
nx = numel(X);

% can we do something more elegant than just min and max? let's instead
% make the range symmetrical around the median, so that we always have a
% single bin (regardless of width) centered around the median 

% also thought about using midhinge instead of median, but median is more
% accepted as a location statisitc, so it makes sense to use that as our
% "baseline"

% first, although we are using histc, let's treat it such that the outer
% edges of the bin always contain min(X) and max(X). We do this via X+1 and X-1

small = 1e-10;
valu = max(X+small)-median(X);
vald = median(X)-min(X-small);

val = max(valu,vald);

minx = median(X) - val;
maxx = median(X) + val; % does this negatively affect how we count series made of integers? well, we don't really have to worry about that at this point

% starting condition only
bins = linspace(minx,maxx,2); 
y = histc(X,bins); 
raty = y / nx; 

for b = 4:2:16 % this will produce b + 1 counts, with the final count indicating events in X that match bins(end)
       
    if max(raty) <= hist_ratio_thr 
        break
        % break, and use the *previous* result (and bins) for the next step
    else    
        % go on to the next iteration; only updates if we have a good result
        nbins = b; % always use an EVEN number of bins so that we always center a bin AROUND the median!
        bins = linspace(minx,maxx,nbins);  
        y = histc(X,bins);
        raty = y / nx;
    end
    
end

% find the limits of this bin
%find(raty == max(raty))
%raty
ind = max(find(raty == max(raty))); % needs to be an actual bin; first find the highest count; if tied, then take the higher value 
% (i.e., this will account for a longer duration; 100 events at 500 ms will take more time than 100 events at 250 ms)
  
minck = bins(ind);
maxck = bins(ind+1); % this will always be valid, since the final bin in the series will always have a count of 0 since it is 1e-10 + max(X)!

xl = X(X>=minck);
xh = X(X<=maxck);

xkeep = intersect(xl,xh);

% now take the median of this bin, and this is the new central tendency
% statistic

location = median(xkeep);

output.min = min(X);
output.mean = mean(X);
output.median = median(X);
output.mode = mode(X);
output.max = max(X);
output.location = location;
output.bins = numel(bins);
