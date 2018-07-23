% compare pairs of trnd distributions and do K-S tests, to make sure that
% we always get null

N = input(' Enter the number of voxels: ');
df = input(' Enter the d.f. of the sample: ');
iter = input(' Enter the number of Monte Carlo iterations: ');

vola = trnd(df,N,1);

ks_res = zeros(iter,2);

for i = 1:iter
    
volb = trnd(df,N,1);

[x, y] = kstest2(vola,volb);

ks_res(i,1) = x;
ks_res(i,2) = y;

end

ks_95 = prctile(ks_res(:,2),5);

ks_sig = sum(ks_res(:,1)) / iter;

figure
hist(ks_res(:,2),10);