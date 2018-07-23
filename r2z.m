function [z] = r2z(r)

% Fisher transform
% note: Fisher z is not the same as the normal distrbution, so
% p-values cannot directly be calculated from this. For the latter,
% see rje code "r2p.m"
%
% rje, july 2012

% Fisher
z = 0.5 * log((1+r) ./(1-r));



