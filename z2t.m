function [tvals] = t2z(zvals,df)

% takes a matrix of t-values and converts them to z-values using MATLAB
% functions

pvals = normcdf(zvals);
tvals = tinv(pvals,df);

