function [Z P] = rosenthal_z(r_ab1, r_ab2, r_b1b2, N, K)

% function [Z P] = rosenthal_z(r_ab1, r_ab2, r_b1b2, N, K)
% formula from Meng et al. 1992
% 
% N = number of subjects
% K = number of regressors (K = 1 for simple correlations)
%
% operates either on single values or three vectors of values
%
% rje, 2012

% 0. get a vector of N if we are using vectors

numv = size(r_ab1,1);

N = zeros(numv,1) + N;
K = zeros(numv,1) + K;

% 1. fisher z-transforms

z_ab1 = 0.5 * log((1+r_ab1) ./(1-r_ab1));
z_ab2 = 0.5 * log((1+r_ab2)./(1-r_ab2));

% 2. f and r2 term

r2 = (power(r_ab1,2) + power(r_ab2,2))/2;

f = (1 - r_b1b2) ./ (2*(1-r2));

f(f > 1) = 1;

% 3. h term

h = (1 - f .* r2) ./ (1 - r2);

% 4. the final equation

Z = (z_ab1 - z_ab2) .* sqrt((N - K - 2) ./ (2 .* (1-r_b1b2) .* h));

% two-tailed p-value

P = (1 - normcdf(abs(Z)))*2;



