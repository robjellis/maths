function output = plot_ecdf(x,y)

% converts x and y to ECDF values and returns correlations

%% exclude NaNs

z = [x y];

% find nan
znan = isnan(z);
znan = max(znan,[],2);

z = z(znan == 0,:);

x = z(:,1);
y = z(:,2);

%% ecdf - values will be in same order as in original input
plot_it = 0;
fx = ecdf_mod(x,plot_it);
fy = ecdf_mod(y,plot_it);

%% plot it

figure(95)
lx = log10(x) / log10(max(x));
ly = log10(y) / log10(max(y));
plot(lx,ly,'.')

figure(99)
% normalized by max 
xnorm = x / max(x);
ynorm = y / max(y);
plot(xnorm,ynorm,'.')

figure(100)
plot(fx,fy,'.')
xlabel('ECDF of x-values')
ylabel('ECDF of y-values')

figure(101)
plot(norminv(fx),norminv(fy),'.')
xlabel('ECDF-to-Z of x-values')
ylabel('ECDF-to-Z of y-values')

%% correlations
spearman = corr(fx,fy,'type','Spearman');
pearson = corr(fx,fy,'type','Pearson');

output.pearson = pearson;
output.spearman = spearman;