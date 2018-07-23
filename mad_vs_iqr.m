function [out] = mad_vs_iqr(num_iti, iter)

% compare MAD vs IQR for potential data

mads = zeros(iter,1);
iqrs = zeros(iter,1);

for i = 1:iter
   
    x = (randn(num_iti,1)+5) * 100;
    
    mads(i) = mad(x) / median(x) * 100;
    iqrs(i) = (iqr(x) / 2) / median(x) * 100; % so that it scales better with mad
    
end

figure(20)

plot(mads,iqrs,'.')


out.corr = corr(mads,iqrs);