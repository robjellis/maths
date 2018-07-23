function output = plot_mahal(x,y,logit,thr)

% converts x and y to ECDF values and returns correlations
% thr (from 0 to 1) is the percentile of observed Mahalnobis distances to use as a
% criteria for outlier detection

%% exclude NaNs

z = [x y];

% find nan
znan = isnan(z);
znan = max(znan,[],2);

z = z(znan == 0,:);

% log transform?
if logit == 1
    z = log10(z);
end

x = z(:,1);
y = z(:,2);



%% calculate mahalnobis
m = mahal(z,z);

% find critical values - all values are >= 0
crit = prctile(m,100*thr);

mout = m >= crit;


%% plot it

figure(200)
plot(x,y,'.')
xlabel('x-values')
ylabel('y-values')

% add in the outliers
xout = x(mout);
yout = y(mout);

hold on
plot(xout,yout,'r.')
hold off


%% correlations
spearman = corr(x,y,'type','Spearman');
pearson = corr(x,y,'type','Pearson');

output.pearson = pearson;
output.spearman = spearman;