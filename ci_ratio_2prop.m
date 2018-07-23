function out = ci_ratio_2prop(a,m,b,n,width)

% formulas are consistent with the literature
%
%
%
%
% RJE | 15 Jan 2018

% note: these scores will always be >= 0 since they *exponentiated in base e*
% to convert them to "signed" values that makes sense in base 2, we simply
% take log2 of the final values


alpha = 1 - width/100;
z = norminv(1-alpha/2); % will be the positive z-value

%% Katz formula
phi = (a/m) / (b/n);

var_hat = (1/a) - (1/m) + (1/b) - (1/n);
lb = exp(log(phi) - z*sqrt(var_hat));
ub = exp(log(phi) + z*sqrt(var_hat));

out.katz_phi = exp(log(phi));
out.katz_ci  = [lb ub];

% alt version of Katz formulas
% RJ confirmed that these yield the same results as above, thus we don't
% need to compute this
do_katz2 = 0;

if do_katz2 == 1
    p1 = a/m;
    q1 = 1 - p1;
    p2 = b/n;
    q2 = 1 - p2;

    phi_katz = p1/p2;

    L = phi_katz * exp(-z*sqrt((q1/(p1*m))+(q2/(p2*n))));
    U = phi_katz * exp(+z*sqrt((q1/(p1*m))+(q2/(p2*n))));

    out.katz2_phi = phi_katz;
    out.katz2_ci  = [L U];
end

%% Adjusted formula

phi_adj = exp(log((a+.5)/(m+.5)) - log((b+.5)/(n+.5)));


% manual corrections only apply to the estimated variance term
if a == m 
    a = a - 0.5;
end

if b == n
    b = b - 0.5;
end

var_hat = (1/(a+.5)) - (1/(m+0.5)) + (1/(b+0.5)) - (1/(n+0.5));

% final
lb = phi_adj * exp(-z*sqrt(var_hat));
ub = phi_adj * exp( z*sqrt(var_hat));

% note that all values below will be >= 0 since we are exponentiating the
% values
out.adj_phi = phi;
out.adj_ci = [lb ub];

% we can also express this in "simpler" terms by taking the log of all
% values; this is more useful when the sign is meaningul
% ACTUALLY, ARE WE ALLOWED TO DO THIS??
out.xxx = '------';

int = [log2(lb) log2(ub)];

% the "score" we care about is the interval endpoint which is closer to
% zero, or zero itself if the interval contains zero

if min(int) > 0
    % if interval is entirely above 0
    score = min(int);
elseif max(int) < 0
    % if interval is entirely below 0
    score = max(int);
else
    % interval contains zero
    score = 0;
end

% final correction - we always want to be careful with very small counts
if a == 0 && b == 0
    score = NaN; % expresss it as NaN so we know that this is intentional
end

out.log2_adj_phi = log2(phi);
out.log2_adj_ci  = int;
out.log2_score   = score;


% or we can express everything in log2 terms
% actually, not sure if this is "allowable"
% phi_log2 = 2.^(log2((a+.5)/(m+.5)) - log2((b+.5)/(n+.5)))
% lb_log2 = phi_log2 * 2.^(-z*sqrt(var_hat))
% ub_log2 = phi_log2 * 2.^( z*sqrt(var_hat))
% 
% phi_alt = (log2((a+.5)/(m+.5)) - log2((b+.5)/(n+.5)))
% lb_alt  = 2.^(phi_alt - z*sqrt(var_hat))
% ub_alt  = 2.^(phi_alt + z*sqrt(var_hat))

%% Bailey formula - looks weird
% p1 = a/m;
% p2 = b/n;
% phi = (a/m) / (b/n);
% 
% num   = z * sqrt(p1/a + p2/b - ((z^2*p1*p2)/(9*a*b))/3)
% denom = 1 - (z^2*p2 / 9*b);
% 
% lb = 
% ub = 