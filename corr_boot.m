function out = corr_boot(x,y,spacing,nbins,subsamp,iter,CI,plot_fig)

% bootstrap correlation coefficient (both Pearson and Spearman) on subsamples from x and y and
% calculate the CI of least squared trend line

if nargin <= 5
    plot_fig = 1;
    CI = 95;
elseif nargin <= 6
    CI = 95;
end
% will also calculate the skewness of x and y, separately

if numel(x) ~= numel(y)
    error('Problem here!')
end

% get rid of NaNs to be safe
nanx = isnan(x);
nany = isnan(y);
nanz = sum([nanx nany],2); % sum the rows

z = nanz >= 1; % if either x or y is NaN

x(z==1,:) = []; % remove these rows
y(z==1,:) = []; % remove these rows

nx = numel(x);

% get x-axis range
xmin = min(x);
xmax = max(x);

plotymin = min(y);
plotymax = max(y);

range = xmax - xmin;
xvals = xmin : range/100 : xmax;
numx = numel(xvals);

r_p    = nan(iter,1);
r_s    = nan(iter,1);

yvals       = nan(iter,numx);
skew_x      = nan(iter,1);
skew_y      = nan(iter,1);

% additional options; get equal values in x-bins?
if strcmp(spacing,'log') || strcmp(spacing,'lin')

    if strcmp(spacing,'lin')
        edges = linspace(xmin,xmax+xmax/1000,nbins);
    elseif strcmp(spacing,'log')
        % note: only words for positive values!
        edges = logspace_rje(xmin,xmax+xmax/1000,nbins); % 100 distinct x-values; make sure xmax is actually larger than observed so the highest value is in the last "real" bin
    end
    
    edges = unique(edges); % don't need to round, but we do only want unique values
    
    nedges = numel(edges);
    
    % figure out what is the smallest number of values in any bin
    hist_res = histc(x,edges); 
    samp = min(hist_res(1:nedges - 1)); % ignore the final bin
    samp = floor(samp / 2); % so we don't use all the values on each permutation
    
    out.xrange = [xmin xmax];
    out.edges = edges;
    out.samp_min = samp;
    
    fprintf([' The minimum number of values in ecah of ' num2str(edges) ' is ' num2str(samp) '.\n'])
end

progressbar(0)

for i = 1:iter

    if strcmp(spacing,'log') || strcmp(spacing,'lin')
        
        xsub = nan(samp,nedges-1);
        ysub = nan(samp,nedges-1);   

        for j = 1:nedges-1
            keep = x >= edges(j) & x < edges(j+1);
            xposs = x(keep);
            yposs = y(keep);
            
            % permute to be safe so we don't take lowest LFIDs
            permind = randperm(sum(keep));
            
            xsub(:,j) = xposs(permind(1:samp));
            ysub(:,j) = yposs(permind(1:samp));
            
        end
        
    else
        % easy
        % random index pairs, with replacement
        inds = randint(subsamp,1,[1 nx]);

        xsub = x(inds);
        ysub = y(inds);       
    end

    % get back into a vector
    xsub = xsub(:);
    ysub = ysub(:);
    
%     if plot_fig == 1
%         figure(10)
%         plot(xsub,ysub,'.')
%     end
    
    % skewness
    skew_x(i) = skewness(xsub,0);
    skew_y(i) = skewness(ysub,0);
    
    % correlations
    r_p(i) = corr(xsub,ysub,'type','Pearson');
    r_s(i) = corr(xsub,ysub,'type','Spearman');
    
    % linear trend
    p = polyfit(xsub,ysub,1);
    
    yvals(i,:) = polyval(p,xvals);
    
    progressbar(i/iter)
end

progressbar(1)

% get CI interval
int_lo = (100 - CI)/2;
int_hi = 100 - int_lo;

ymax = prctile(yvals,int_hi);
ymed = prctile(yvals,50);
ymin = prctile(yvals,int_lo);

% get CI of Pearson and Spearman
pearson_CI  = prctile(r_p, [int_lo int_hi]);
spearman_CI = prctile(r_s, [int_lo int_hi]);

% CI for skewness
skew_x_CI = prctile(skew_x,[int_lo int_hi]);
skew_y_CI = prctile(skew_y,[int_lo int_hi]);

%% plot

if plot_fig > 0
    figure(100)

    % full data
    full_r_p = corr(x,y,'type','Pearson');
    full_r_s = corr(x,y,'type','Spearman');
    % plot the original data
    plot(x,y,'.')
    xlabel('x')
    ylabel('y')
    hold on
    
    % show the bins
    if nbins > 2
        for j = 1:nedges
            plot([edges(j) edges(j)],[plotymin plotymax],'m')
        end
    end
    
%     % the trend for the full data 
%     p = polyfit(x,y,1);
%     y = polyval(p,xvals);

    plot(xvals,ymed,'r','LineWidth',1.5)
    plot(xvals,ymax,'r:')
    plot(xvals,ymin,'r:')

    hold off

    % in fact, the ECDF plots just fill up the screen and make it hard to
    % see anything!
    
%     % now plot with ECDF values; we have already cleared out NaNs
     figure(101)
     plot(x,ecdf_mod(y),'.')
    
    % skewness results
    figure(102)
    subplot(1,2,1)
    ecdf(skew_x)

    subplot(1,2,2)
    ecdf(skew_y)
    
    

    
end

%% output

out.num_pairs = numel(x);
out.full_r_p = full_r_p;
out.full_r_s = full_r_s;
out.pearson_CI = pearson_CI;
out.spearman_CI = spearman_CI;
out.skew_x_CI = skew_x_CI;
out.skew_y_CI = skew_y_CI;
