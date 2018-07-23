function output = corr_perm(x,y,iter)

% permutation test on correlation
%

%{

Other resources:
http://www.mathworks.com/matlabcentral/fileexchange/34920-correlation-permutation-test-with-correction-for-multiple-comparisons/content/mult_comp_perm_corr.m
http://www.uvm.edu/~dhowell/StatPages/Resampling/RandomCorr/randomization_Correlation.html
%}



x = x(:);
y = y(:);

% make sure data is same size
if numel(x) == numel(y)   
    % OK
    n = numel(x);
else
    return
end

% set up variables
corrP = nan(iter,1);
corrS = nan(iter,1); 

% loop
progressbar(0)

for i = 1:iter
    % randomly permute y and keep x the same
    yind = randperm(n);
    yperm = y(yind); % resorts it
    
    % get Pearson and Spearman correlations
    corrP(i) = corr(x,yperm,'type','Pearson');
    corrS(i) = corr(x,yperm,'type','Spearman');
    
    progressbar(i/iter)
end

progressbar(1)

% now get the distribution
[pf px] = ecdf_mod(corrP);
[sf sx] = ecdf_mod(corrS);

figure(100)
plot(x,y,'.')

figure(101)
subplot(1,2,1)
plot(px,pf,'.')
ylabel('Pearson r')

subplot(1,2,2)
plot(sx,sf,'.')
ylabel('Spearman rho')

% now return the p-value of the *actual* correlation

[rP pP] = corr(x,y,'type','Pearson');
[rS pS] = corr(x,y,'type','Spearman');

% abs value
arP = abs(rP);
arS = abs(rS);

% find where the observed r-value falls on the ecdf
ecdf_P = min(pf(px >= rP));
ecdf_S = min(sf(sx >= rS));


%% outputs
output.corrP_obs = rP;
output.ecdf_P = ecdf_P;

output.corrS_obs = rS;
output.ecdf_S = ecdf_S;
