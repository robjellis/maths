function mc_sd(N,iter)

% 
% mc_sd(N,iter)
% 
% question: how many events do we need to get a reasonable estimate of the
% SD of a population of values with a known SD [randn]

% use a normal distribution; should be reasonable approximation for gait

sds = zeros(iter,N); % this will have the SD "trace" (columns) for each iteration 

for i = 1:iter
vals = randn(N,1); 

% since we are using N(0,1), we know that the target SD is in fact 1.0
% how close do we get to that?

for j = 1:(N-1)
    sds(i,j) = std(vals(1:j+1));
end

end

figure(300)
for i=1:iter
   plot(sds(i,1:N-1))  
   hold on
end


% now plot the 95% percentiles and mean
plot(prctile_nist(sds(:,1:N-1),5),'r','LineWidth',2);
plot(prctile_nist(sds(:,1:N-1),95),'r','LineWidth',2);
plot(prctile_nist(sds(:,1:N-1),50),'g','LineWidth',2); % median

hold off
