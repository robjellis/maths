function ccdf(x,mult)

% RJE version of this function because Matlab seemed to be slightly different
% "mult" allows percentage rather than proportion to be returned

if nargin == 1
    mult = 1;
end

x = x(:);

n = numel(x);

% sort it
x = sort(x);

% the ccdf values
ccdf = nan(n,1);

for i = 1:n
    ccdf(i) = mult * sum(x > x(i)) / n;
end

figure(100)
plot(x,ccdf)
xlabel('x')
ylabel('1 - CDF')

output.x = x;
output.ccdf = ccdf;