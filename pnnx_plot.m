function pnnx_plot(x,inc)

% RJE version of this function because Matlab seemed to be slightly different
% "mult" allows percentage rather than proportion to be returned

if nargin == 1
    inc = 1;

end
  
x = x(:);

% for PNN values, we just care about 0 to 100

y = 0:inc:100;

n = numel(y);

% sort it
x = sort(x);

% the ccdf values
ccdf = nan(n,1);

for i = 1:n
    ccdf(i) = 100 * sum(x > y(i)) / numel(x);
end

figure(100)
hold on
plot(y,ccdf)

% plot pNN20
%hold on
%plot(y(21),ccdf(21),'.')
%hold off

xlabel('x')
ylabel('pNNx')

output.y = x;
output.ccdf = ccdf;