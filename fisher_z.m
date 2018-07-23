function [z p] = fisher_z(r1,r2,n1,n2)

% function [z p] = fisher_z(r1,r2,n1,n2)
% Fisher z test
%
% see, e.g., http://www.vassarstats.net/rdiff.html

% first convert r to z
z1 = 0.5 * log((1+r1) ./(1-r1));
z2 = 0.5 * log((1+r2) ./(1-r2));

% now the z-test
z = (z1 - z2) / (sqrt( 1/(n1-3) + 1/(n2-3)));

% now a 2-tailed p-value (divide by 2 for 1-tailed)

p = 2*normcdf(-abs(z));