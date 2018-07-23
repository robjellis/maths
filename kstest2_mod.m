function output = kstest2_mod(samp1,samp2)

% This is a modified function by RJE that will perform a *two-tailed* KS
% test and return an accurate Z-value (using two-tailed p-value) that
% relfects the direction of the test
%
% A positive Z means samp1 has a higher distribution of values than samp2