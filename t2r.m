function [r df] = t2r(t,nsub,nreg)

% [r df] = t2r(t,nsub,nreg)
%
% Standard formula for converting r-value to t-value.
% For a simple correlation, nreg = 1.
%
% df is output just as confirmation that it is as expected
%
% rje, june 2012


r = t ./ sqrt(nsub - nreg - 1 + t.^2);

df = nsub - nreg - 1;