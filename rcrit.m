function [rcrit] = rcrit(nsub,nreg, alpha, tail)

% critical value for Pearson's r
%
% [rcrit] = rcrit(nsub,nreg, alpha, tail)
%
% nsub  = number of subjects
% nreg  = number of regressors (1 for simple correlation)
% alpha = the desired p-value (can be a vector)
% tail  = 1 for one-tail (positive), 2 for two-tail
%
% by rje

% checked against Altman (1991) Table B7

df = nsub - nreg - 1;

if tail == 1
   pcrit = 1 - alpha;
elseif tail == 2
   pcrit = 1 - alpha/2; 
end

% formula:
% http://en.wikipedia.org/wiki/Pearson_product-moment_correlation_coefficient#Testing_using_Student.27s_t-distribution

tcrit = tinv(pcrit,df);
rcrit = tcrit ./ sqrt(df + tcrit.^2);
