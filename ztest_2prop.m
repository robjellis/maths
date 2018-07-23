function out = ztest_2prop(x1,n1,x2,n2)

% RJE | Jan 2018

% compute the "raw" z-value, don't worry about p-values for now
% will also work when x1 and x2 are [N x 1] vectors of counts, and n1 and n2 are single values

% note: if both x1 and x2 == 0, then z will be NaN; this allows for later
% identification and handling (e.g., convert to zero) if desired.

p1 = x1 ./ n1;
p2 = x2 ./ n2;

p  = (x1 + x2) ./ (n1 + n2);

z  = (p1 - p2) ./ sqrt(p.*(1 - p) .* (1./n1 + 1./n2));


% modified version by RJE; keep the observed ratio of x2/n2, but adjust it to match
% sample size of n1; the resultant z_mod value will be more conservative (i.e., closer to zero) 

x3 = (x2 ./ n2) .* n1; % it's ok if not a whole number
n3 = n1;

p3 = x3 ./ n3;

p  = (x1 + x3) ./ (n1 + n3);

z_mod  = (p1 - p3) ./ sqrt(p.*(1 - p) .* (1./n1 + 1./n3));


out.z     = z;
out.z_mod = z_mod;
