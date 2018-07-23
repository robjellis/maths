% bootstrapping to compare observed r-values between a fixed x-variable
% (IV) and random permutations of a y-variable

% x1, x2, x3, x4 are the regressors: Age, TD, MD-perf, RD-perf,
% respectively, put into an X-by-4 array

iter = input(' How many iterations to perform?: ');

numy = numel(x(:,1));

rvals = zeros(iter,4);

for i = 1:iter
   y = randn(numy,1);
for j = 1:4

   rvals(i,j) = corr(x(:,j),y(:,1));
end

end

x1_x2 = corr(rvals(:,1),rvals(:,2))
x2_x3 = corr(rvals(:,2),rvals(:,3))
x3_x4 = corr(rvals(:,3),rvals(:,4))