function [z] = dunn_clark_z(r_ab1, r_ab2, r_b1b2, N)

% formula taken from Hittner 2003
% This is the actual D and C (1969) formula (not the modification by Steiger), 
% as described in Hittner
%
% assumes a single collumn of each value!


% 0. get a vector of N if we are using vectors

numv = size(r_ab1,1);

N = zeros(numv,1) + N;

% 1. fisher z-transforms

z_ab1 = 0.5 * log((1+r_ab1) ./(1-r_ab1));
z_ab2 = 0.5 * log((1+r_ab2)./(1-r_ab2));

% 2. get the cov term (call it "X")

X = ( r_b1b2 .* (1 - r_ab1.^2 - r_ab2.^2) - 0.5 * (r_ab1.*r_ab2) .* (1 - r_ab1.^2 - r_ab2.^2 - r_b1b2.^2) ) / ((1 - r_ab1.^2).*(1 - r_ab2.^2));

% 3. the equation
z = sqrt(N - 3).*(z_ab1 - z_ab2) .* (2 - 2 * X).^(-1/2);