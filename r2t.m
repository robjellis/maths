function [t] = r2t(r,nsub,nreg)

% [t] = r2t(r,nsub,nreg)
%
% standard formula for converting r-value to t-value
% for a simple correlation, nreg = 1
%
% rje, jun 2012



t = r .* sqrt((nsub - nreg - 1) ./ (1- r.^2));