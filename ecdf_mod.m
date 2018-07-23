function [f x] = ecdf_mod(data, plot_it,color)

% ECDF, but returning values in the same order they appear in "data", to
% make transforms easier. This is accomplished using simple resorting of
% matrix values.
%
% * modified in Feb 2016 to properly handle NaNs
% 
% output will contain two fields:
%  x is the original values in data, in the original order
%  f is the corresponding ECDF value


if nargin < 2
    plot_it = 0;
end

if nargin < 3
    color = 'b';
end

% just to be safe
data = data(:);

% sort by value
[data_sort orig_ind] = sort(data); % NaN will be at the end

% find NaN indices
nan_inds = find(isnan(data_sort) == 1);
keep_inds = find(isnan(data_sort) == 0);

% keep only actual values
data_true = data_sort(keep_inds);

[e x] = ecdf(data_true);

% x values will be unique(data), but have a double value of the lowest x

e = e(2:end);
x = x(2:end);

%% figure
if plot_it == 1
    figure(5)
    plot(x,e,color)
    xlabel('x')
    ylabel('Cumulative probability')
    title(' Empirical CDF')
end

% now we need to match up data with x and sub in e

% first, we do unique on the original true data
[data_unique xxx ind2] = unique(data_true); % orig_ind is the key way we get back to the original order

% confirm this
if sum(data_unique - x) == 0
    % perfect
else
    error('Problem with re-sorting of values')
end

% get back to sorted order
x_resort = data_unique(ind2);
f_resort = e(ind2);
    
% get the NaNs back
x_resort(nan_inds) = NaN;
f_resort(nan_inds) = NaN;

% and restore the original order; this has to be done carefully

Z = [orig_ind x_resort f_resort];

Z = sortrows(Z,1); % sort by the original index

x = Z(:,2);
f = Z(:,3);


