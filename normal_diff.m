function [output] = normal_diff(voxels,iter)

dim = round(voxels^(1/3));
vol_mean = zeros(iter,1);
vol_var = zeros(iter,1);
ks_pval = zeros(iter,1);

for i = 1:iter
    vol = randn(dim,dim,dim);
    vol_flip = flipdim(vol,1);

    vol_diff = vol - vol_flip;
    vol_diff = vol_diff(:);
    
    % mean and variance
    vol_mean(i) = mean(vol_diff);
    vol_var(i) = var(vol_diff);

    % don't do kstest, because it will compare against *standard* normal
    % [i.e., N(0,1)]
    
    % kstest
    %[h p] = kstest(vol_diff);
    %ks_pval(i) = p;
    
end

figure(20)
hist(vol_var,30)

output.vol_means = median(vol_mean);
output.vol_vars = median(vol_var);
%output.ks_alpha = sum(ks_pval < .05)/iter;
