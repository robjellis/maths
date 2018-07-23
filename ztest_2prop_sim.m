function out = ztest_2prop_sim(nA,nB,weight,iter)

% nice way of showing how, when there is a greater agreement between
% relative counts in A and B, the observed z-distribution approaches normality

% first get the correlation between counts in A and B
x = randn(iter,1);
y = x*weight + randn(iter,1)*(1-weight);

% center
x = x - mean(x);
y = y - mean(y);

% convert to CDF values
xp = normcdf(x);
yp = normcdf(y);

% get the random counts

cA = round(xp * nA);
cB = round(yp * nB);

% cA = randi(nA,iter,1);
% cB = randi(nB,iter,1);

% z-test on proportions
res = ztest_2prop(cA,nA,cB,nB);

z = res.z;

% ignore NaNs
z_nan = sum(isnan(z));
z(isnan(z)) = [];

figure(19)
plot(cA/nA,cB/nB,'.')

figure(20)
ecdf(z)

figure(21)
qqplot(z)

%% output
out.rho     = corr(cA,cB,'type','Spearman');
out.z_mn    = mean(z);
out.z_std   = std(z);
out.z_skew  = skewness(z,0);
out.z_nan   = z_nan;
