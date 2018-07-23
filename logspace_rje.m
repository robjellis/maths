function out = logspace_rje(a,b,n)

%
% logspace_rje(a,b,n)
%
% log spacing between any values a and b, in n steps

loga = log(a);
logb = log(b);

% now use linspace
exps = linspace(loga,logb,n);

% restore
out = exp(exps);
