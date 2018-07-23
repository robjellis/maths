function [stats] = ba_calc_new(vec1,vec2,plot_it,fignum,data1_name,data2_name,ba_mod)

% Bland-Altman plots and stats
%
% vec1 and vec2 are vectors of identical shape; if either has an NaN case,
% that case will be automatically DROPPED
%
% see http://en.wikipedia.org/wiki/Bland%E2%80%93Altman_plot
%
% Since the BA calculation calculates mean deviation as "vec1 - vec2", vec1
% should be the "experimental/novel", and vec2 should be the "control"
%
% RJE | 2013.03.21
% updated with regression model 2014.07.09
%

if nargin<3
   plot_it = 1;
end

if nargin < 4
   fignum = 20;
end

if nargin < 6
   data1_name = 'data1';
   data2_name = 'data2';
end

if nargin < 7
    ba_mod = 0; % ba_mod == 1 means slightly modified x-axis: just use vec2
end
    
vec1 = vec1(:);
vec2 = vec2(:);

%% need to check if there are NaN, and remove them
% revised: if either x or y has an NaN value, that point will not be
% plotted; thus we don't need to exclude NaN cases; furthermore, by
% excluding these cases, we MESS UP the indexing of the original input
% series!

c1 = isnan(vec1);
c2 = isnan(vec2);
csum = c1+c2;

vec1_tmp = vec1(csum == 0); % use this for CALCULATING, not PLOTTING
vec2_tmp = vec2(csum == 0); % use this for CALCULATING, not PLOTTING

if numel(vec1) ~= numel(vec2)
   % bad
   fprintf('\n Error: vec1 and vec2 do not have the same shape.\n')
   return
else
   % ok
end

% for plotting
if ba_mod == 0
    baX = (vec1 + vec2) / 2;
elseif ba_mod == 1
    baX = vec2;
end

baY =  vec1 - vec2; % this will affect where the MEAN value falls (either positive or negative)
                    % rje thinks it makes sense to do 
                    % "novel measure - standard measure"

%% standard model                    
% for calculating
baX_tmp = (vec1_tmp + vec2_tmp) / 2;
baY_tmp =  vec1_tmp - vec2_tmp; 

ba_mean = mean(baY_tmp);
ba_2sd  = 2 * std(baY_tmp); % use 2 rather than 1.96 per Bland Altman 1999

% below is method by RJE
% * additional thought: why not use KDE to estimate the shape of the PDF,
% and then get the limits that way!
ba_lo =  prctile_nist(baY_tmp, 2.3); % equivalent to - 2SD
ba_md =  prctile_nist(baY_tmp, 50); % median
ba_hi = prctile_nist(baY_tmp, 97.7); % equivalent to + 2SD

% how many y-axis points are > 2SD?
ba_mean_tmp = ba_mean(isnan(ba_mean) == 0); % manually ignore nans
abs_err = abs(baY_tmp - ba_mean_tmp);

error_rate = sum(abs_err > ba_2sd) / numel(vec1_tmp) * 100; % error rate percentage


%% regression mode

% Bland Altman 1999 just use simple linear regression
    % step 1: get the coefficients
    coef = polyfit(baX_tmp,baY_tmp,1);
    
    % step 2: get regression line
    y_est = (coef(1) * baX_tmp) + coef(2);

    % step 3: get the residuals
    resid = baY_tmp - y_est;
    resid_mn = mean(resid);
    resid_sd = std(resid);
    
    % step 4: now, we calculate the upper and lower lines
    resid_hi = y_est + 2 * resid_sd;
    resid_lo = y_est - 2 * resid_sd;
    
    % for sake of comparison, let's fit KDE to this too
    npoints = 2^12;
    min_res = min(resid);
    max_res = max(resid);
    ran_res = max(resid) - min(resid);
    
    minv = min_res - ran_res/4;
    maxv = max_res + ran_res/4;
    
    % we take the smaller bw estimate of Botev KDE and Matlab ksdensity
    [bot_bw, xxx, xmesh] = kde(resid,npoints,minv,maxv);
    [xxx, xxx, mat_bw] = ksdensity(resid,xmesh,'npoints',npoints);
    
    min_bw1 = min(bot_bw,mat_bw);
    
    % use ksdensity just to be sure things fit
    pdf = ksdensity(resid,xmesh,'npoints',npoints,'width',min_bw1);
    
    % make sure min = 0
    pdf = pdf - min(pdf);
    
    % scale to make max = 1.0 (arbitrary)
    pdf = pdf / max(pdf);
    
    % now get the cdf (still in arbitrary units)
    cdf = cumsum(pdf);
    
    % scale to 1.0
    cdf = cdf / max(cdf);
       
    % originally tried to directly get 2SD, but if N is small then tails may be poorly estimated.
    % instead, we will get the estimated 1 SD, and then multiply by 2
    kde_med = max(xmesh(cdf <= 0.5));
    kde_1n  = max(xmesh(cdf <= .1587));
    kde_1p  = min(xmesh(cdf >= .8413)); 
    
    % do this carefully so we automatically handle distributions with positive values
    kde_lo  = kde_med + 2 * (kde_1n - kde_med);
    kde_hi  = kde_med + 2 * (kde_1p - kde_med);
    
    % now we need to get the plottable lines
    reg_kde_hi = kde_hi - resid_mn;
    reg_kde_lo = resid_mn - kde_lo;  
   

% RJE: why not use robust regression in combination with kernel density estimation!

    % 1. use robust regression
       [rob_coef stats] = robustfit(baX_tmp,baY_tmp);
       p = stats.p;
       rob_slope_p = p(2);
       
       % get the best_fit line
       rob_est = (rob_coef(2) * baX_tmp) + rob_coef(1);
       
       % get the residuals
       rob_res = baY_tmp - rob_est;
       rob_mn = mean(rob_res);
       rob_sd = std(rob_res); % we DON'T use this because we will use KDE

    % 2. use Botev KDE to fit the data
    min_rob = min(rob_res);
    max_rob = max(rob_res);
    ran_rob = max(rob_res) - min(rob_res);
    
    minv = min_rob - ran_rob/2;
    maxv = max_rob + ran_rob/2;
    
    [bot_bw, xxx, xmesh] = kde(rob_res,npoints,minv,maxv);
    [xxx, xxx, mat_bw]   = ksdensity(rob_res,xmesh,'npoints',npoints);
    
    min_bw2 = min(bot_bw,mat_bw);
    
    % use ksdensity just to be sure things fit
    pdf = ksdensity(rob_res,xmesh,'npoints',npoints,'width',min_bw2);
    
    % make sure min = 0
    pdf = pdf - min(pdf);
    
    % scale to make max = 1.0 (arbitrary)
    pdf = pdf / max(pdf);
    
    % now get the cdf (still in arbitrary units)
    cdf = cumsum(pdf);
    
    % scale to 1.0
    cdf = cdf / max(cdf);
       
    % originally tried to directly get 2SD, but if N is small then tails may be poorly estimated.
    % instead, we will get the estimated 1 SD, and then multiply by 2
    
    rob_1sd_lo = max(xmesh(cdf <= .1587));
    rob_1sd_hi = min(xmesh(cdf >= .8413));
    
    rob_lo = rob_1sd_lo * 2;
    rob_hi = rob_1sd_hi * 2;
    
    % now we need to get the plottable lines
    plus_err = rob_hi - rob_mn;
    min_err  = rob_mn - rob_lo;
  
    figure(22)
    plot(xmesh,pdf)
    hold on
    plot(xmesh,cdf,'r')
    
    % now plot the 1 SD value
    plot([rob_1sd_lo rob_1sd_lo],[0 1],'g')
    plot([rob_1sd_hi rob_1sd_hi],[0 1],'g')
    hold off



%% plot

if plot_it == 1
    
    if fignum > 0
       figure(fignum)
    else
        % don't call a new figure
    end


    %plot(baX,baY,'*') % "*" tends to work better since "+" can blend in with axis ticks
                      % we use the ORIGINAL vector because this is the only
                      % way that NaN cases will be preserved
                      
    plot(baX,baY,'Marker','.','MarkerSize',10,'LineStyle','none')
    hold on
    
    % BA error bars
    min_y = ba_mean - ba_2sd;
    max_y = ba_mean + ba_2sd;

    %min_x = floor(min(baX));
    %max_x =  ceil(max(baX));
    min_x = min(baX);
    max_x = max(baX);

    plot([min_x max_x],[min_y min_y],'b');
    plot([min_x max_x],[ba_mean ba_mean],'b','LineStyle','--');
    plot([min_x max_x],[max_y max_y],'b');

    % percentiles
    %plot([min_x max_x],[ba_lo ba_lo],'m');
    %plot([min_x max_x],[ba_md ba_md],'m','LineStyle','--');
    %plot([min_x max_x],[ba_hi ba_hi],'m');    
    
    % plot the standard regression results
    plot(baX_tmp, y_est,'r','LineWidth',1,'LineStyle','--')
    plot(baX_tmp, resid_hi,'r','LineWidth',1)
    plot(baX_tmp, resid_lo,'r','LineWidth',1)
    
    % plot the KDE results for comparison
    plot(baX_tmp, y_est - reg_kde_lo,'r','LineWidth',1,'LineStyle','--')
    plot(baX_tmp, y_est + reg_kde_hi,'r','LineWidth',1,'LineStyle','--')
    
    % plot robust regression results
   
    plot(baX_tmp, rob_est,'m','LineWidth',1,'LineStyle','--')
    plot(baX_tmp, rob_est + plus_err,'m','LineWidth',1)
    plot(baX_tmp, rob_est - min_err,'m','LineWidth',1)    
    
    hold off

    if ba_mod == 0
        lx = xlabel(['(' data1_name ' + ' data2_name ') / 2']);
    elseif ba_mod == 1
        lx = xlabel(data2_name);
    end   
    ly = ylabel([data1_name ' - ' data2_name]);
    title(['Bland-Altman plot' sprintf('\n') '2*SD = ' num2str(ba_2sd) '; > 2*SD = ' num2str(error_rate) '%']);
    set(lx,'Interpreter','none'); set(ly,'Interpreter','none');
    
    % residual plots
        
    figure(21)
    % regular residuals
    plot(baX_tmp,resid,'r+')
    hold on
    
    % robust residuals
    plot(baX_tmp,rob_res,'mo')
    
    % regular error bars
    lo_val = resid_mn - 2 * resid_sd;
    hi_val = resid_mn + 2 * resid_sd;
    plot([min_x max_x], [lo_val lo_val], 'r');
    plot([min_x max_x], [resid_mn resid_mn],'r','LineStyle','--')
    plot([min_x max_x], [hi_val hi_val], 'r');
    
    plot([min_x max_x], [kde_lo kde_lo], 'r','LineStyle','--')
    plot([min_x max_x], [kde_hi kde_hi], 'r','LineStyle','--');   
    
    % robust + KDE error bars
    plot([min_x max_x], [rob_lo rob_lo], 'm');
    plot([min_x max_x], [rob_mn rob_mn], 'm','LineStyle','--');
    plot([min_x max_x], [rob_hi rob_hi], 'm');
    
    hold off


end

% spearman correlation on BA plot
[ba_rho ba_pval] = corr(baX,baY,'type','Spearman');    
 %% outputs
stats.data_pearson  = corr(vec1,vec2);
stats.data_spearman = corr(vec1,vec2,'type','Spearman');
stats.BA_mean  = ba_mean;
stats.BA_2sd = ba_2sd;
stats.BA_rho = ba_rho;
stats.BA_pval = ba_pval;
stats.error_rate = error_rate;
stats.BA_perc_lo = ba_lo;
stats.BA_perc_md = ba_md;
stats.BA_perc_hi = ba_hi;
stats.BA_perc_lo_diff = ba_md - ba_lo;
stats.BA_perc_hi_diff = ba_hi - ba_md;