function out = bino_ci_calc(x,n,width)
%
% bino_ci_calc(x,n,width)
%
% Compute popular CIs for binomial proportions of x successes out of n trials.
% width is the desired width/confidence; e.g., 80%, 90%, 95%
%
% Data is entered as actual counts, because that's how the formulas take them
% note: fractional numbers are allowed, and the math works fine on them!
%
% all methods tested have a manual correction for x == 0 or x == n
%
% References consulted:
% * Vollset 1993 - the easiest formulas
% * Newcombe 1998
% * Agresti & Coull 1998 - their approximation is  not computed here, as
%   it's only valid for a 95% CI
%
% RJE | 2017.11.18

% only works on single values for now
if numel(x) > 1 || numel(n) > 1
    error(' Only single values of x and n are allowed.')
end

if nargin < 3
    width = 95;
end

% initialize
wald_nc = [];
wald_cc = [];
wils_nc = [];
wils_cc = [];
exact   = [];

% key other variables
alpha = 1 - width/100;
z = norminv(1- alpha/2);

p = x / n;

% compute the CI - be careful with the ( ) ! Newcombe paper actually
% writes the formulas in the same way that they can be typed in Matlab

%% A. Simple or Wald method (no continuity correction)
% Newcome 1998 formulas are equivalent to Vollset
% wald_nc_lo = p - z*sqrt(p*q/n);
% wald_nc_hi = p + z*sqrt(p*q/n);

if x == 0 % adjustment, per Vollset
    wald_nc_lo = 0;
    wald_nc_hi = 1 - exp(log(alpha/2)/n);

elseif x == n % adjustment, per Vollset
    wald_nc_lo = exp(log(alpha/2)/n);
    wald_nc_hi = 1;

else
    % Formulas from Vollset 1993
    wald_nc_lo = (x/n) - z/sqrt(n)*sqrt((x/n)*(1-x/n));
    wald_nc_hi = (x/n) + z/sqrt(n)*sqrt((x/n)*(1-x/n));
end

% additional logical correction by RJE
if wald_nc_lo < 0
    wald_nc_lo = 0;
elseif wald_nc_hi > 1
    wald_nc_hi = 1;
end

wald_nc = [wald_nc_lo wald_nc_hi];



%% B. Simple or Wald method (with continuity correction)
% Formula from: Newcombe 1998 - equivalent to Vollset
% wald_cc_lo = p - (z*sqrt(p*q/n) + 1/(2*n));
% wald_cc_hi = p + (z*sqrt(p*q/n) + 1/(2*n)) ; 

if x == 0 % adjustment, per Vollset
    wald_cc_lo = 0;
    wald_cc_hi = 1 - exp(log(alpha/2)/n);

elseif x == n % adjustment, per Vollset
    wald_cc_lo = exp(log(alpha/2)/n);
    wald_cc_hi = 1;

else % Formulas from Vollset 1993
    wald_cc_lo = (x/n) - (z/sqrt(n) * sqrt((x/n)*(1- x/n)) + 1/(2*n));
    wald_cc_hi = (x/n) + (z/sqrt(n) * sqrt((x/n)*(1- x/n)) + 1/(2*n));         
end

% additional logical correction by RJE
if wald_cc_lo < 0
    wald_cc_lo = 0;
elseif wald_cc_hi > 1
    wald_cc_hi = 1;
end

wald_cc = [wald_cc_lo wald_cc_hi];


%% C. Wilson score, no continuity correction
% Based on Agresti & Coull eq. 2
% wilson_nc_lo = (p + z^2/(2*n) - z*sqrt((p*q + z^2/(4*n))/n)) / (1 + z^2/n);
% wilson_nc_hi = (p + z^2/(2*n) + z*sqrt((p*q + z^2/(4*n))/n)) / (1 + z^2/n);

% Formulas from: Vollset 1993 (uses x rather than p and q)
% Matches results from R package PropCI

% RJ note: the special if/elseif/else setup prevents some slight hiccups in
% Matlab (e.g., a value of wilson_nc_hi = 1.0000000000)

if x == 0
    wilson_nc_lo = 0; % logical correction
    wilson_nc_hi = (x + z^2/2 + z*sqrt(x - x^2/n + z^2/4)) / (n + z^2);
    
elseif x == n 
    wilson_nc_lo = (x + z^2/2 - z*sqrt(x - x^2/n + z^2/4)) / (n + z^2);
    wilson_nc_hi = 1; % logical correction
    
else % standard case
    wilson_nc_lo = (x + z^2/2 - z*sqrt(x - x^2/n + z^2/4)) / (n + z^2);
    wilson_nc_hi = (x + z^2/2 + z*sqrt(x - x^2/n + z^2/4)) / (n + z^2);
end

wils_nc = [wilson_nc_lo wilson_nc_hi];

%% D. Wilson score, with continuity correction

% Formula from: Vollset 1993; simply add or subtract 0.5 from the success count (x)
wilson_cc_lo = ((x - 0.5) + z^2/2 - z*sqrt((x - 0.5) - (x - 0.5)^2/n + z^2/4)) / (n + z^2);
wilson_cc_hi = ((x + 0.5) + z^2/2 + z*sqrt((x + 0.5) - (x + 0.5)^2/n + z^2/4)) / (n + z^2);

% adjustment per Newcombe
if p == 0
    wilson_cc_lo = 0;
elseif p == 1
    wilson_cc_hi = 1;
end

wils_cc = [wilson_cc_lo wilson_cc_hi];



%% E. Clopper-Pearson (Exact) - don't worry about this (also slows down the code)
% Formulas from Agresti 1998 (also gives adjustments for x == 0 and x == n)
% Baker 2000 notes that this formula comes from Lemmis & Trivedi 1996
% see also: http://jansenlex.readyhosting.com/mwsug/2008/pharma/MWSUG-2008-P08.pdf for clear formulas
% formulas from Leemis 1996 Appendix C don't seem to work

% % Formulas from  Agresti 1998 - doesn't seem to work
% Flo = finv((1-alpha)/2 , 2*x    , 2*(n-x+1));
% Fhi = finv(alpha/2    , 2*(x+1), 2*(n-x));
% 
% ex_lo = 1 / (1 + (n-x+1) /      x*Flo);
% ex_hi = 1 / (1 + (n - x) / ((x+1)*Fhi));
% 
% exact = [ex_lo ex_hi];
% 
% % correction for boundaries
% if x == 0
% elseif x == n
% end

%% outputs
out.interval = width;
out.observed_p = x/n;
out.wald_nc = wald_nc;
out.wald_cc = wald_cc;
out.wils_nc = wils_nc;
out.wils_cc = wils_cc;
out.exact   = exact;
