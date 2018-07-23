function [stats] = plot_eda(data1,data2,opt,fignum,data1_name,data2_name,ba_mod,do_plot)
%
% plot_eda(data1,data2,opt,fignum,figname)
%
% Note: data1 should be "standard" and data2 should be "experimental
%
% a single function that allows user to plot any or all of 1D and 2D simple
% plot. data2 can be set to [] and will be ignored, yielding only 1D plots
%
% "opt" can be a single number or set in [ ]; e.g., [3 4 5]:
%
% 1 ECDF (1 or 2 D)
% 2 histogram (1 or 2 D)
% 3 Box plot (1 or 2 D)
% 4 Scatter plot (2D only)
% 5 Bland-Altman plot (2D only)
% 6 Quantile-Quantile plot (2D only)
%
% data1 and data2 can be any dimension, but will be sorted to a single column;
% thus, they must have the same dimensional structure (and handedness);
% this is important for 3D neuroimaging data
%
% "figname" allows a helpful label to be applied so figures are easy to
% tell apart
% 
% see also:
% http://www.mathworks.com/products/statistics/description3.html
% http://www.mathworks.com/help/stats/scatterhist.html


if nargin < 2 % decide if 1 or 2 sets of dat
   ndata = 1;
   data2 = [];
   opt = [1 2 3]; % these are the only things we can do
   
elseif numel(data2) == 0 % keep this here, in case we have the input plot_eda(x,[],[1 2 3]) etc
   ndata = 1;
   opt = [1 2 3];    
else
   ndata = 2;
end

if nargin < 3 && ndata == 2
   opt = [1 2 3 4 5 6];
end

nopt = numel(opt);

if max(opt) > 6
   opt = [1 2 3 4 5 6];
end

if nargin < 4
   fignum = 10;
end

if nargin < 5
   data1_name = 'Data1';
   data2_name = 'Data2';
end

if nargin < 7
    ba_mod = 0; % "1" means just use data2 as x-axis
    do_plot = 1;
end

if nopt == 1
    r = 1;
    c = 1;
elseif nopt == 2
    r = 1;
    c = 2;
elseif nopt < 4
    r = 1;
    c = 3;
elseif nopt <=6 
    r = 2;
    c = 3;
end


%% data size stuff

if ndata == 1
    data1 = data1(:);
    data1_tmp = data1(isnan(data1) == 0);
    nx_tmp = numel(data1_tmp);
    
    % get min and max
    minV = min(data1);
    maxV = max(data1);

elseif ndata == 2
    data1 = data1(:);
    data2 = data2(:);
    
    % pull out NaN
    keep_ind = and(isnan(data1) == 0,isnan(data2) == 0);
    data1 = data1(keep_ind);
    data2 = data2(keep_ind);
    
    if sum(size(data1) - size(data2)) == 0
        % OK
    else
        fprintf(' Error: data1 and data2 are not the same size.\n\n')
    end

    % If either x or y has an NaN value, that point will not be
    % plotted; thus we don't need to exclude NaN cases; furthermore, by
    % excluding these cases, we MESS UP the indexing of the original input
    % series!

    % Revised: we want to exclude cases with NaN for sorted plots and box
    % plots, but *not* Scatter plots or B-A plots (since this will mess up case
    % numbers)

    c1 = isnan(data1);
    c2 = isnan(data2);
    csum = c1 + c2;

    data1_tmp = data1(csum == 0);
    data2_tmp = data2(csum == 0);

    nx_tmp = numel(data1_tmp);
    
    % get min and max
    min1 = min(data1);
    max1 = max(data1);
    min2 = min(data2);
    max2 = max(data2);
    
    % get global min and max
    minV = min(min1,min2);
    maxV = max(max1,max2);
    
end

if do_plot == 1
    figure1 = figure(fignum); % take the number from param
    
    clf % will clear all previous content from this figure
        
    if ismember(1,opt) || ismember(2,opt) || ismember(3,opt)
        % Create textbox
        annotation(figure1,'textbox',[0.012 0.97 0.06781 0.03681],...
            'String',{data1_name},...
            'Interpreter','none',...
            'FontWeight','bold',...
            'FontSize',11,...
            'FitBoxToText','off',...
            'EdgeColor','none',...
            'Color',[1 0 0]);

        if ndata == 2
            % Create textbox
            annotation(figure1,'textbox',[0.012 0.93 0.06781 0.03681],...
                'String',{data2_name},...
                'Interpreter','none',...
                'FontWeight','bold',...
                'FontSize',11,...
                'FitBoxToText','off',...
                'EdgeColor','none',...
                'Color',[0 0 1]);
        end
    end
end


%% Work
for p = 1:nopt
        subplot(r,c,p)
        choice = opt(p);
        plot_this(choice) % subfunction
end

%output = 'Done';
    

%% subfunction: plot_this
    
    function plot_this(choice)
    
        if choice == 1 % ECDFs; exclude NaNs
            
            if ndata == 1 || ndata == 2
               [f1 x1] = ecdf_rje(data1_tmp); % ecdf_rje is non-proprietary
               stairs(x1,f1,'color','r','LineWidth',2)
            end
            
            if ndata == 2
                [f2 x2] = ecdf_rje(data2_tmp);
                hold on
                stairs(x2,f2,'color','b','LineWidth',1)
                hold off
            end
            
            lx = xlabel('Value'); ly = ylabel('Cumulative probability'); 
            set(lx,'Interpreter','none'); set(ly,'Interpreter','none');
            title('ECDF');
            
            % axis limits
            axis([minV maxV 0 1])
            
        
        elseif choice == 2 % Botev's KDE estimation
%             % set number of bins using Rice's rule
%             % https://en.wikipedia.org/wiki/Histogram#Rice_Rule
%             nbins = ceil(2 * nx_tmp^(1/3));
%             
%             ncases = round(nx_tmp / nbins); % number of cases per bin
%             
%             if ndata == 1
%                x = linspace(min(data1_tmp),max(data1_tmp),nbins);
%             elseif ndata == 2
%                x = linspace(min(min(data1_tmp),min(data2_tmp)),max(max(data1_tmp),max(data2_tmp)),nbins);
%             end
%             
%             
%             if ndata == 1
%                y1 = histc(data1_tmp,x);
%                plot(x,y1,'LineWidth',1,'color','r','DisplayName',data1_name,'LineWidth',2)
%             
%             
%             elseif ndata == 2 
%                y1 = histc(data1_tmp,x);
%                y2 = histc(data2_tmp,x);
%                
%                plot(x,y1,'LineWidth',1,'color','r','DisplayName',data1_name,'LineWidth',2)
%                hold on
%                plot(x,y2,'LineWidth',1,'color','b','DisplayName',data2_name)
%                
%                hold off
%        
%             end

            if ndata == 1
                [b density xmesh] = kde(data1_tmp); % use defaults
                plot(xmesh,density,'LineWidth',1,'color','r','DisplayName',data1_name,'LineWidth',2)
            elseif ndata == 2
                % need to set a common range for the two plots (pulled from kde.m)
                Range = maxV - minV;
                MIN = minV-Range/10; 
                MAX = maxV+Range/10; 
                
                [b density xmesh] = kde(data1_tmp,2^12,MIN,MAX);
                plot(xmesh,density,'LineWidth',1,'color','r','DisplayName',data1_name,'LineWidth',2)
                hold on
                
                [b density xmesh] = kde(data2_tmp,2^12,MIN,MAX);
                plot(xmesh,density,'LineWidth',1,'color','b','DisplayName',data2_name,'LineWidth',1)
                hold off
            end
            
            %lx = xlabel('Value'); ly = ylabel('Count'); 
            %title(['Histogram outlines' sprintf('\n(') num2str(ncases) ' cases per bin)'])
            
            lx = xlabel('Value');
            ly = ylabel('Density');
            
            set(lx,'Interpreter','none'); set(ly,'Interpreter','none');
            
            title('Kernel density estimation')
        
        elseif choice == 3 % box plot
               if ndata == 1
                   data = data1;
               elseif ndata == 2
                   data = [data1, data2]; % columns
               end

               percentiles = [1 99];
               boxplot_rje(data,percentiles,0); % rje function (free, no stats toolbox required)

               hold off
       
            lx = xlabel('Data'); ly = ylabel('Value'); 
            set(lx,'Interpreter','none'); set(ly,'Interpreter','none');
            title('Box plot')
            % note: axis are set in boxplot_rje

            
        elseif choice == 4 % scatter plot; preserve NaN cases
            plot(data1,data2,'*')
            lx = xlabel(data1_name); ly = ylabel(data2_name); 
            set(lx,'Interpreter','none'); set(ly,'Interpreter','none');
            
            % manually calculate the Pearson r-value (so we don't need "corr.m" from the stats toolbox
            pearson = cov(data1,data2) / (std(data1)*std(data2)); % per http://en.wikipedia.org/wiki/Pearson_product-moment_correlation_coefficient
            pearson = pearson(1,2);

            
            % now plot the reference line (code borrowed from Matlab qqplot)

            q1x = prctile_nist(data1,25);
            q2x = prctile_nist(data1,50);
            q3x = prctile_nist(data1,75);
            
            q1y = prctile_nist(data2,25);
            q2y = prctile_nist(data2,50);
            q3y = prctile_nist(data2,75);
%             
%             qx = [q1x; q3x];
%             qy = [q1y; q3y];
% 
%             dx = q3x - q1x;
%             dy = q3y - q1y;
%             slope = dy./dx;
%             centerx = (q1x + q3x)/2;
%             centery = (q1y + q3y)/2;
%             
%             maxx = prctile_nist(data1,.05); % avoid extreme outliers
%             minx = prctile_nist(data1,99.5);
%             
%             maxy = centery + slope.*(maxx - centerx); % passes through center
%             miny = centery - slope.*(centerx - minx);
% 
%             mx = [minx; maxx];
%             my = [miny; maxy];

            % get the regression
            coef = polyfit(data1,data2,1);

            % step 2: get regression line
            mx = [min1 max1];
            
            y1 = coef(1) * min1 + coef(2);
            y2 = coef(1) * max1 + coef(2);
            
            my = [y1 y2];
            
            hold on
            plot(mx,my,'Color','m','LineStyle','-','Marker','none','LineWidth',1.5);
            hold off
            
            title(['Scatter plot' sprintf('\n') 'Pearson r = ' num2str(pearson)]);
            
            % axis limits
            axis([min1 max1 min2 max2]) % keep the box "tight" so it's easier to see the data
            
        elseif choice == 5 % Bland-Altman plot; preserve NaN cases
            plot_ba = 1;
            fignum = 0;
            [ba_stats] = ba_calc(data2,data1,plot_ba,fignum,data2_name,data1_name,ba_mod); 
            
            % axis limits
            xx = (data2+data1)/2;
            yy = data2-data1;
            ymin = min(yy);
            ymax = max(yy);
            xmin = min(xx);
            xmax = max(xx);
            
            % in case all y-axis values are identical (i.e., 0) ...
            if ymin == ymax
               ymin = ymin - 1;
               ymax = ymax + 1;
            end
            
            axis([xmin xmax ymin ymax])
            
        elseif choice == 6 % Q-Q plot; ignore NaN cases
            qqplot_rje(data1_tmp,data2_tmp,50,0,data1_name,data2_name); % "50" is the number of quantiles; "0" means don't call a new figure
            
        end

    end
set(gcf, 'Color', 'w');

%% data cursor

if do_plot == 1
    datacursormode on;
    h = datacursormode(figure1);
    set(h,'UpdateFcn',@myupdatefcn);
end

    function [indtxt] = myupdatefcn(obj,event_obj)
        pos = get(event_obj,'Position');

        index = get(event_obj,'DataIndex');

        indtxt = ['Ind: ' num2str(index) sprintf('\n') '  X: ' num2str(pos(1)) sprintf('\n') '  Y: ' num2str(pos(2))];
    end

%% stats output
stats.data1_mean = mean(data1);
stats.data1_sd = std(data1);

if nopt >= 4
    stats.data1_p25 = q1x;
    stats.data1_p50 = q2x;
    stats.data1_p75 = q3x;
end

if numel(data2) > 0
    stats.break1 = 'xxxxxx';
    stats.data2_mean = mean(data2);
    stats.data2_sd = std(data2);
    stats.break2 = 'xxxxxx';
    stats.BA_mean = ba_stats.BA_mean;
    stats.BA_2SD = ba_stats.BA_2sd;
    stats.BA_rho = ba_stats.BA_rho;
    stats.BA_perc_lo = ba_stats.BA_perc_lo;
    stats.BA_perc_md = ba_stats.BA_perc_md;
    stats.BA_perc_hi = ba_stats.BA_perc_hi;
    stats.BA_perc_lo_diff = ba_stats.BA_perc_lo_diff;
    stats.BA_perc_hi_diff = ba_stats.BA_perc_hi_diff;    
    
    if nopt >=4
        stats.data2_p25 = q1y;
        stats.data2_p50 = q2y;
        stats.data2_p75 = q3y;

    end
end

end