function s = skewness_vec(x,bias)

% RJE function for sample skewness
% * requires a vector ([1 x N] or [N x 1])
% * will scrub all NaNs out before calculating
%
% based on Matlab function skewness.m
 
%   S = SKEWNESS(X) returns the sample skewness of the values in X.  For a
%   vector input, S is the third central moment of X, divided by the cube
%   of its standard deviation.  For a matrix input, S is a row vector
%   containing the sample skewness of each column of X.  For N-D arrays,
%   SKEWNESS operates along the first non-singleton dimension.
%
%   SKEWNESS(X,0) adjusts the skewness for bias.  SKEWNESS(X,1) is the same
%   as SKEWNESS(X), and does not adjust for bias.
%

% confirm that we have a vector
if size(x,1) == 1 || size(x,2) == 1
    % OK
else
    % quit
    fprintf('\n Error: skewness_vec only works on vector format.\n\n')
    return
end

if nargin < 2
   bias = 1;
end

x = x(isnan(x) == 0); % get rid of NaN

% Center x, compute its third and second moments, and compute the
% uncorrected skewness.
x0 = x - mean(x);
s2 = mean(x0.^2); % this is the biased variance estimator
m3 = mean(x0.^3);
s = m3 ./ s2.^(1.5);

% Bias correct the skewness.
if bias == 0
    n = sum(~isnan(x));
    n(n<3) = NaN; % bias correction is not defined for n < 3.
    s = s .* sqrt((n-1)./n) .* n./(n-2);
end
