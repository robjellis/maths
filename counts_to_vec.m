function V = counts_to_vec(C,T)

% turn a vector of counts into a vector of those counts
% C = counts, in order; e.g., C = [20 10 10 40 60]
% T = target values, in order; e.g., T = [1 2 3 4 5]
%
% if nargin = 1, assume T = 1:numel(C)
%
% RJE 11 Dec 2017

if nargin == 1
    T = 1:numel(C);
end

% how many unique values?
nV = numel(T);

% how many total items?
tot = sum(C);

V = nan(tot,1);
ctr = 1;

for i = 1:nV
    
    V(ctr:ctr+C(i)-1) = repmat(T(i),C(i),1);
    
    ctr = ctr + C(i);
    
end

