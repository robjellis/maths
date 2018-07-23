% qq plot simulation

N = input(' Number of voxels: ');
df = input(' Enter the d.f.: ');
iter = input(' Enter the number of iterations: ');
%steps = input(' Enter the number of quantiles: ');

t = -12:12;

%steps = 0:2:100;
steps = 100 * (1 ./ (1 + exp(-1*t)));
%steps = [.0001 .0005 .001 .005 .01 .05 1 5 10 25 50 75 90 95 99 99.95 99.99 99.995 99.999 99.9995 99.9999];
prt = input(' Enter the percentile to evaluate: ');

%z = linspace(0,100,steps);
qcorr = zeros(iter,1);

for i = 1:iter
   vol1 = trnd(df,N,1);
   vol2 = trnd(df,N,1);
   
   qx = prctile(vol1,z)';
   qy = prctile(vol2,z)';
   
   qcorr(i) = corr(qx,qy);
    
end

%figure
%hist(qcorr,30);

figure
qqplot(vol1,vol2,steps);

prt_val = prctile(qcorr,prt)
