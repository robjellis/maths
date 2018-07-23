function [stats] = ba_calc(vec1,vec2,plot_it,fignum,data1_name,data2_name,ba_mod)

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
    ba_mod = 0; % slightly modified x-axis: just use vec2
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

vec1_tmp = vec1(csum == 0); % use this for CALCULATING (can't have NaN), not PLOTTING
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

                    
% for calculating
if ba_mod == 1
    baX_tmp = vec2_tmp;
else
    baX_tmp = (vec1_tmp + vec2_tmp) / 2;
end
baY_tmp =  vec1_tmp - vec2_tmp; 

ba_mean = mean(baY_tmp);
ba_2sd  = 2 * std(baY_tmp); % adjusted to 2 (rather than 1.96) to be a bit more conservative

% below is method by RJE
ba_lo =  prctile_nist(baY_tmp, 2.3); % equivalent to - 2SD
ba_md =  prctile_nist(baY_tmp, 50); % median
ba_hi = prctile_nist(baY_tmp, 97.7); % equivalent to + 2SD

% how many y-axis points are > 2SD?
ba_mean_tmp = ba_mean(isnan(ba_mean) == 0); % manually ignore nans
abs_err = abs(baY_tmp - ba_mean_tmp);

error_rate = sum(abs_err > ba_2sd) / numel(vec1_tmp) * 100; % error rate percentage

%% now calculate the regression line

coef = polyfit(baX_tmp,baY_tmp,1);

% step 2: get regression line
y_est = (coef(1) * baX_tmp) + coef(2);
    
%% plots
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
    
    plot([min_x max_x],[min_y min_y],'r');
    plot([min_x max_x],[ba_mean ba_mean],'r','LineStyle','--');
    plot([min_x max_x],[max_y max_y],'r');

    % percentiles
    %plot([min_x max_x],[ba_lo ba_lo],'m');
    %plot([min_x max_x],[ba_md ba_md],'m','LineStyle','--');
    %plot([min_x max_x],[ba_hi ba_hi],'m');    
    
    % regression line
    plot(baX_tmp,y_est,'m')
    hold off

    if ba_mod == 0
        lx = xlabel(['(' data1_name ' + ' data2_name ') / 2']);
    elseif ba_mod == 1
        lx = xlabel(data2_name);
    end   
    ly = ylabel([data1_name ' - ' data2_name]);
    title(['Bland-Altman plot' sprintf('\n') '2*SD = ' num2str(ba_2sd) '; > 2*SD = ' num2str(error_rate) '%']);
    set(lx,'Interpreter','none'); set(ly,'Interpreter','none');
end
    
 %% outputs
stats.Pearson  = corr(vec1,vec2);
stats.Spearman = corr(vec1,vec2,'type','Spearman');
stats.BA_mean  = ba_mean;
stats.BA_2sd = ba_2sd;
[rho p] = corr(baX,baY,'type','Spearman');
stats.BA_rho = rho;
stats.BA_p = p;
stats.error_rate = error_rate;
stats.BA_perc_lo = ba_lo;
stats.BA_perc_md = ba_md;
stats.BA_perc_hi = ba_hi;
stats.BA_perc_lo_diff = ba_md - ba_lo;
stats.BA_perc_hi_diff = ba_hi - ba_md;