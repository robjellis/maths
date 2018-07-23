function [y] = prctile_nist(data, p)
%
% Compute one or more percentiles from data
% FUNCTION [y] = spm_percentile(data, p)
% data - [N x M]. Each column (M) will have it's own percentile(s) calculated
%        (after excluding NaNs will be excluded)
% p    - scalar or vector of percentage values (from 0 to 100)
%        if not specified, p defaults to all quartiles: [0 25 50 75 100]
%
% y    - scalar or n-vector of corresponding percentiles
%
%__________________________________________________________________________
% Copyright (C) 2010 Wellcome Trust Centre for Neuroimaging
%
% *** Original code by Ged Ridgway (for the SPM Group) ("spm_percentile.m")
% *** modified by Robert J Ellis (http://robjellis.net)
%
% version = 2013.03.06
%

% The algorithm used is that described by NIST, with x = k+d = 1+p(N-1)/100
% http://www.itl.nist.gov/div898/handbook/prc/section2/prc262.htm
%
% This choice apparently matches Excel, but not MATLAB's prctile, though 
% the differences are typically small (and are zero for min, median, max).
%
% This algorithm was chosen because it requires no special handling of 0 or
% 100, and lets 0 and 1 percentiles differ even with less than 100 samples.
% It also has the appealing property of returning uniformly spaced values
% for uniformly spaced percentiles of uniformly spaced data. For example:
%  x = 1:11;
%  p = 0:25:100;
%  diff([prctile(x, p(:)) spm_percentile(x, p) spm_percentile_nist(x, p)])
% gives constant differences only for spm_percentile.
%

data_full = data;

if nargin < 2
    p = 0:25:100;
end
if any(p < 0 | p > 100)
    fprintf('Percentage values outside 0 <= p <= 100 are not allowed')
    return
end

p    = p(:); % will write out as a column
np = numel(p);

nc = size(data,2);

y = nan(np,nc);

for i = 1:nc
    data = data_full(:,i);

    data        = sort(data(~isnan(data))); % get rid of NaN
    N           = length(data);
    
    if N > 0
        data(end+1) = data(end); % to handle special case of 100th percentile below

        x = 1 + p * (N - 1) / 100;
        k = floor(x);
        d = x - k;

        y(:,i) = (1 - d) .* data(k) + d .* data(k + 1);
    else
        % no non-NaN values!
        y(:,i) = NaN;
    end
end
