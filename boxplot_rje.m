function output = boxplot_rje(data, percentiles, fignum,plot_outliers,box_type)

%
% modified by RJE (2013.04.08) from
% http://www.mathworks.com/matlabcentral/fileexchange/15053-simple-box-whiskers-plot-for-those-who-dont-have-the-statistics-toolbox
%
% RJE explictly makes the whiskers be a percentile value: [low high], e.g., percentiles = [1 99]
%
% box_type is a boxplot by default, or 'e' for a simple median-plus-outlier threshold
%
% example:
% x = randn(1000,1);
% boxplot_rje(x,[1 99],0)
%
% updated in Jan 2018 with colors


    %data = sort(data, 1); % ascend
    
    if nargin < 2
       percentiles = [2.5 97.5]; % to yield a 95% interval
    end
    
    if nargin < 3
       fignum = 500;
    end
    
    if nargin < 4
       plot_outliers = 0;
    end
    
    if nargin < 5
        box_type = 'b';
    end
    
    % other defaults
    width = 1;
    lineWidth = .5;


    [m n] = size(data);

    % get the full min and max of the data
    minn = min(min(data));
    maxx = max(max(data));
      
    
    q2 = prctile_nist(data,50); % median
    
    q1 = prctile_nist(data,25); % 25th percentile
    q3 = prctile_nist(data,75); % 75th percentile
    
    min_whisk = prctile_nist(data,percentiles(1)); 
    max_whisk = prctile_nist(data,percentiles(2)); 
    

    if fignum == 0
        % don't make a new figure
        plot(0,0) % just to clear the window
    else

        figure(fignum)  
        plot(0,0) % just to clear the window
    end
    hold on 

    if strcmp(box_type,'b')
        draw_data = [max_whisk; q3; q2; q1; min_whisk];
            
        n = size(draw_data, 2);

        unit = (1-1/(1+n))/(1+9/(width+3));
     
        for i = 1:n

            % calculate outliers
            xx = data(:,i);

            xx1 = xx(xx < min_whisk(i));
            xx2 = xx(xx > max_whisk(i));

            xx = [xx1; xx2;]; % just the outliers, in a column
            nx = numel(xx);

            v = draw_data(:,i);

            % draw the low whisker
            plot([i-unit/2, i+unit/2], [v(5), v(5)], 'LineWidth', lineWidth,'Color',[0 0 0]);

            % draw the high whisker
            plot([i-unit/2, i+unit/2], [v(1), v(1)], 'LineWidth', lineWidth,'Color',[0 0 0]);

            % draw median
            plot([i-unit, i+unit], [v(3), v(3)],'r', 'LineWidth', 1.6,'Color',[1 0 0]); 

            % draw vertical lines connecting whiskers to box
            plot([i, i], [v(5), v(4)], 'LineWidth', lineWidth,'LineStyle','--','Color',[0 0 0]);
            plot([i, i], [v(2), v(1)], 'LineWidth', lineWidth,'LineStyle','--','Color',[0 0 0]);

            % draw box
            plot([i-unit, i+unit, i+unit, i-unit, i-unit], [v(2), v(2), v(4), v(4), v(2)], 'LineWidth', lineWidth,'LineStyle','-','Color',[0 0 0]);
            
            % highlight first and third quartile
            plot([i-unit, i+unit], [v(2), v(2)], 'LineWidth', 1.2,'Color',[0 0 1]);
            plot([i-unit, i+unit], [v(4), v(4)], 'LineWidth', 1.2,'Color',[0 0 1]);

            % plot the outliers
            if plot_outliers == 1
               plot(repmat(i,[nx 1]),xx,'r+')
            end

            % set the axis
            if minn == maxx
                % special weird case
                axis([0 n+1 minn-.01 maxx+.01])
            else
                axis([0 n+1 .99*minn 1.01*maxx])
            end

        end
        
    else
        % simpler plot
        % have to fix the error bars a bit: "When they are vectors, each
        % error bar is a distance of L(i) below and U(i) above the point
        % defined by (X(i),Y(i))." - http://www.mathworks.com/help/matlab/ref/errorbar.html
        L = q2 - min_whisk; % now it is the length
        U = max_whisk - q2; 
        errorbar(1:n,q2,L,U,'.')
    end

    hold off
    set(gcf,'color','w');
%% export stats: percentiles    
out_prc = [2.5 25 50 75 97.5]; 

if size(data,1) > 1
    output = prctile_nist(data,out_prc); % percentiles will be in columns
else
    output = [];
end
    