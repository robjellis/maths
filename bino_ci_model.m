function out = bino_ci_model(p,N,iter)

% Plot ECDF of average success rate when we take random samples of size N from a binomial
% distribution with success rate p


if nargin < 2
    N = power(2,3:6);
end

if nargin < 3
    iter = 5000;
end

% pre-allocate
means = nan(iter,1); % for a single N-level

figure(50)
clf
hold on

for j = 1:numel(N)
   
   for i = 1:iter
       
       vals = rand(N(j),1) <= p; % will yield a vector in which there are 100*p% 1s
       
       means(i) = mean(vals);
       
   end
   
   % add the ECDF
   [f, x] = ecdf(means);
   stairs(x,f)
   
      
end

hold off


    
    