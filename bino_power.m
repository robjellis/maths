function bino_power(p1,p2,alpha,power)

% bino_power(p0,p1,alpha,power)
%
% Compute sample size (for an A/B test) to achieve desired power at 
% target alpha level (two-tailed).
%
% Note that
%   p1 = x and p2 = x + y, versus
%   p1 = x and p2 = x - y
%
% will have *different* resultant sample sizes; this is to be expected.
%
% Sample size without and with continuity correction is calculated.
%
% Based on: Fleiss JL, Tytun A, Ury HK (1980): 
%           A simple approximation for calculating sample sizes for comparing independent proportions. 
%           Biometrics 36:343–6.
%
% See also: http://math.furman.edu/~dcs/courses/math47/R/library/Hmisc/html/bpower.html 
%
% RJE | Nov 2017

if nargin < 3
    alpha = .05;
end

if nargin < 4
    power = .80;
end

d = abs(p2 - p1); % this term is squared in the equation, so sign doesn't matter

q1 = 1 - p1;
q2 = 1 - p2;

P = (p1 + p2)/2;
Q = 1 - P;

z_a = norminv(1 - alpha/2); % must be alpha/2
z_b = norminv(power); % power = 1 - beta

N = (z_a*sqrt(2*P*Q) + z_b*sqrt(p1*q1 + p2*q2))^2 / d^2;

% round up
N = ceil(N);


% continuity correction (optional)
N_cc = (N/4)*(1 + sqrt(1 + 4/(N*d)))^2;

% round up
N_cc = ceil(N_cc);


% results
fprintf(['\n  Each study arm (A and B) must have at least: ',...
         '\n   ' num2str(N) ' subjects (without continuity correction), or',...
         '\n   ' num2str(N_cc) ' subjects (with continuity correction). \n\n']);