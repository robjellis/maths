function output = nancorr(x,y,plot_it)

% to get rid of NaNs, -Inf, and +Inf, in x and y and then do the correlation
%
% by RJE | last updated 2018.01.06

if nargin < 3
    plot_it = 0; % don't plot
end

x = x(:);
y = y(:);

%% get rid of NaNs, +Inf, and -Inf to be safe
nanx = isnan(x);
nany = isnan(y);

negIx = x == -Inf;
negIy = y == -Inf;

posIx = x == Inf;
posIy = y == Inf;

mat = [nanx nany negIx negIy posIx posIy];

keep = max(mat,[],2); % take max value for each row

x = x(keep == 0);
y = y(keep == 0);



%% now do the correlation
[coef1, pval1] = corr(x,y,'type','Pearson');
[coef2, pval2] = corr(x,y,'type','Spearman');

if plot_it == 1
   figure(10)
   subplot(1,2,1)
   plot(x,y,'.') % visualize Pearson
   xlabel('x')
   ylabel('y')
   
   % ranks (visualize Spearman)
   xr = tiedrank(x);
   yr = tiedrank(y);
   subplot(1,2,2)
   plot(xr,yr,'.')
   xlabel('x ranks')
   ylabel('y ranks')
    
end


%% output
output.pearson_coef = coef1;
output.pearson_pval = pval1;
output.spearman_coef = coef2;
output.spearman_pval = pval2;
output.num_pairs = numel(x);




