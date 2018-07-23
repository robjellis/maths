function [output] = corr_simp(N,iter)

r = zeros(N,1);

for i = 1:iter
    
    A = randn(N,1);
    B = randn(N,1);

    r(i) = corr(A,B);

end

figure(10)

hist(r,30)

output.prc95 = prctile(r,95);
output.prc99 = prctile(r,99);
output.prc995 = prctile(r,99.5);


