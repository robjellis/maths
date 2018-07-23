nc = corr_norm(v1,v2);

% the normalized correlation [cost function] of Jenkinson et al. (2002)
%
% rje, august 2012

nc = sum(v1.*v2) / (sqrt(sum(v1.^2))*sqrt(sum(v2.^2)));