function out = outlier_rje(data,method,thr,plot_it)

%
% out = outlier_rje(data,method,thr,plot_it)
%
% method of detecting outliers by RJE
%
% option 1: using central tendency / location ('l')
%   1. find the central tendency; median by default
%   2. get absolute deviations from central tendency
%   3. sort these values
%   4. take the "first-order ratio" of the sorted values
%   5. use "thr" to find the first value that has a ratio > thr
%   6. exclude all *subsequent* values in the sorted series
%   7. return the indices to the user
% option 2: using median absolute deiviation ('mad')
%   1. based on the assumption that STD is ~ 1.5 x MAD for normal distribution
%   2. Find the MAD
%   3. Multiply this value x 1.4826 *and* the user-specified threshold to get
%      the equivalent "number of STDs above the median" (default = 3)
% option 3: using standard deviation ('sd')
%   1. This is the simplest method, but may be biased by outliers
%   2. Defined outliers as falling outside the range [mean - thr * SD, mean + thr * SD]

mult = 1.5;
% mult = 1.4826; % technically, this is only valid for large samples (> 5000)

% defaults

if nargin < 2
    method = 'mad';
end

if nargin < 3
    thr = 3;
end

if nargin < 4
    plot_it = 0;
end

data = data(:); % just to be safe

% 0. are we calculating using TEMPORALLY ORDERED data or not? 
if strcmp(method,'l')
    % do nothing; data stays as is
    
elseif strcmp(method,'s') % we identify outliers based on temporal order
    % we take first order difference of raw data
    data = diff(data);
    
    % restore the index
    data = [0; data];
end

% 2. get the abs diffs
if strcmp(method,'l') || strcmp(method,'mad')
    loc = nanmedian(data); % need to do this in case there are cases with NaN that we need to preserve for sake of temporal order etc.
    diffs = abs(data - loc);
elseif strcmp(method,'sd')
    loc = nanmean(data);
    diffs = abs(data - loc);
elseif strcmp(method,'s')
    diffs = abs(data);
end

ind = 1:numel(diffs);

% 3. sort these, and make sure we stay paired with the ORIGINAL DATA

data(:,2) = ind;
data(:,3) = diffs;
sorted = sortrows(data,3);

data_sort = sorted(:,1);
ind_sort = sorted(:,2);
diffs_sort = sorted(:,3);

% 4a. first-order ratio
rat = diffs_sort(2:end) ./ diffs_sort(1:end-1);

% 4b. adjust the index location
rat = [0 rat'];


% 5. find the outliers

if strcmp(method,'mad')
    
    % need a little protection here in case median is very small ...
    diffs_mad = prctile_nist(diffs_sort,50);
    diffs_mad = max(diffs_mad,.003);
    rat(diffs_sort < diffs_mad) = 0;
    
    thr_final = thr * mult * diffs_mad;
    ind_outliers = ind(diffs_sort > thr_final); 
    % mult to convert MAD to STD; 
    % 3 is to look for values > 3x larger than STD, based on
    % http://en.wikipedia.org/wiki/68%E2%80%9395%E2%80%9399.7_rule
elseif strcmp(method,'sd')
    thr_final = thr * nanstd(data(:,1));
    ind_outliers = ind(diffs_sort > thr_final);
else 
    ind_outliers = ind(rat > thr); % this will just be the actual indices, not 0s and 1s
end

% 6. exclude *subsequent* indices too
ind_first = min(ind_outliers);
ind_outliers = ind_first:numel(diffs_sort);

% 7. isolate the original values and indices
ind_out = ind_sort(ind_outliers);
val_out = data_sort(ind_outliers);

if plot_it == 1
    
    figure(305)
    clf
    
    subplot(3,1,1)
    plot(ind,data(:,1),'.')    
    
    subplot(3,1,2)
    plot(ind,diffs_sort,'+')
    
    if strcmp(method,'mad')
        % show the MAD and the final threshold
        hold on
        plot([1 numel(diffs_sort)],[diffs_mad diffs_mad],'g')
        hold on
        plot([1 numel(diffs_sort)],[thr_final thr_final],'r')
        hold off
    end

    subplot(3,1,3)
    if strcmp(method,'mad')
        % leave blank
    else
        plot(1:numel(diffs_sort),rat,'m+')
        hold on
        plot([1 numel(diffs_sort)],[thr thr],'g')
        hold off
    end

    % plot on the data
    subplot(3,1,1)
    hold on
    plot(ind_out,val_out,'ro')
    hold off
    
    title(['Method: ' method])
    
    % back to first subplot
    subplot(3,1,2)
    hold on
    plot(ind(ind_outliers),diffs_sort(ind_outliers),'ro')
    hold off
    
end

%% output

out.inds = ind_out; % the indices of the outliers
out.vals = val_out; % the values of the outliers
out.thr_lo = loc - thr_final;
out.thr_hi = loc + thr_final;
out.perc = 100 * numel(ind_out) / numel(data(:,1)); % percentage of input data that has been flagged as an outlier