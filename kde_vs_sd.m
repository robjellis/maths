function output = kde_vs_sd(n,iter)

x=randn(n,1);

% simple estimate
mn = mean(x);
sd = std(x);


mn_1n = mn - sd
mn_1p = mn + sd

% kde estimate

% for sake of comparison, let's fit KDE to this too
npoints = 2^12;
min_x = min(x);
max_x = max(x);
ran_x = max(x) - min(x);

minv = min_x - ran_x;
maxv = max_x + ran_x;

minv = -6;
maxv = 6;

% we take the smaller bw estimate of Botev KDE and Matlab ksdensity
[bot_bw, xxx, xmesh] = kde(x,npoints,minv,maxv);
[xxx, xxx, mat_bw]   = ksdensity(x,xmesh,'npoints',npoints);

min_bw = min(bot_bw,mat_bw);

% use ksdensity just to be sure things fit
pdf = ksdensity(x,xmesh,'npoints',npoints,'width',min_bw);

% make sure min = 0
pdf = pdf - min(pdf);

% scale to make max = 1.0 (arbitrary)
pdf = pdf / max(pdf);

% now get the cdf (still in arbitrary units)
cdf = cumsum(pdf);

% scale to 1.0
cdf = cdf / max(cdf);

figure(10)
plot(x,ones(size(x))*.1,'.')
hold on
plot(xmesh,pdf)
hold on
plot(xmesh,cdf,'r')
hold off

% originally tried to directly get 2SD, but if N is small then tails may be poorly estimated.
% instead, we will get the estimated 1 SD, and then multiply by 2
kde_med = max(xmesh(cdf <= 0.5));

% 1 SD
kde_1n  = max(xmesh(cdf <= .1587)) 
kde_1p  = min(xmesh(cdf >= .8413)) 

% 1.96 SD to yield 95%
kde_2n  = max(xmesh(cdf <= .025));
kde_2p  = min(xmesh(cdf >= .975));

% do this carefully so we automatically handle distributions with positive values
kde_lo  = kde_med + 2 * (kde_1n - kde_med);
kde_hi  = kde_med + 2 * (kde_1p - kde_med);

