function [p t z1 z2] = r2p(rho,tail,nsub,nreg)

% function [p t z1 z2] = r2p(rho,tail,nsub,nreg)
%
% returns the p-value for Pearson r-value
% code copied from corr.m in MATLAB ("pvalPearson"), and modified slightly by rje
% MATLAB: t(k) = rho(k).*sqrt((n-2)./(1-rho(k).^2));
%
% z1 is actual r to z (normal distribution)
% z2 is Fisher's z
%
% tail is 'b' 'r' or 'l'
% nsub is number of subjects
% nreg is number of regressors (will be 1 for a simple correlation)
%
% these p-values agree with Statistica results, and with Altman (1991)
% Table B7 (note: he only considers a partial correlation with one dummy
% variable, and thus suggests an n (- 1) - 1 correction
%
% rje, july 2012

%t = sign(rho) .* Inf
%k = (abs(rho) < 1)
t = rho .* sqrt((nsub-nreg-1)./(1 - (rho .^2)));

% t2z for good measure
p = tcdf(t,nsub-nreg-1);
z1 = norminv(p);

% Fisher z, for comparison
z2 = 0.5 * log((1+rho) ./(1-rho));

switch tail
case 'b' % 'both or 'ne'
    p = 2*tcdf(-abs(t),nsub-nreg-1);
case 'r' % 'right' or 'gt'
    p = tcdf(-t,nsub-nreg-1);
case 'l' % 'left' or 'lt'
    p = tcdf(t,nsub-nreg-1);
end