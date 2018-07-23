function output = robustfit_slope(x,y,split_method,param,fit_method,plot_it,fig_num)

% calculate the slope, as a percentage change, using a window with size pts
% which *preceeds* the value (to make it more like actual perception)
%
% 'split_method' is either 'all' (individual points) or 'sec' (for time-based sections)
% 'fit_method' is either use robustfit ('r') *or* simple polyfit ('p')
%
% RJE

if nargin < 5
    plot_it = 0;
end

if nargin < 6
    fig_num = 200;
end

% *************************
% optional shift by half a window
shift = 1;
% *************************
        
% size of data
x = x(:);
y = y(:);

% initialize slope_pc
slope_pc = NaN;

nvals = numel(y);

% split_method choice
if strcmp(split_method,'all')
    % slope based on "param" preceeding points, at *all* points
    slope_pc = nan(nvals,1);
    pts = param;
    for i = pts:nvals

        yv = y(i-pts+1:i);

        if isempty(x)
            xv = cumsum(yv);
        else
            xv = x(i-pts+1:i);
        end

        b = robustfit(xv,yv);

       % find the change along the Y axis relative to intercept
       s2 = b(2)*xv(end) + b(1); % b(1) is the intercept value
       s1 = b(2)*xv(1)   + b(1);

       % 2013.09.04 - percentage change from y1 to y2 (the current point)
       slope_pc(i) = 100 * (s2 - s1) / s1; 

    end
elseif strcmp(split_method,'sec')
    % slope calculated for "param"-s sections
    
    slice = param;
     
        if isempty(x)
            x = cumsum(y);
        else
            % we use existing y
        end
    
        minx = min(x);
        maxx = max(x);

        % original position
        sections = minx:slice:maxx; 
        sections = sections(:);
        starts = sections(1:end-1);
        ends = sections(2:end); 

        
        if shift == 1
            % shift by half a window
            sections2 = (minx+slice/2):slice:maxx; 
            sections2 = sections2(:);
            starts2 = sections2(1:end-1);
            ends2 = sections2(2:end);  
        else
            starts2 = [];
            ends2 = [];
        end
        
        % merge them
        starts = [starts; starts2];
        ends = [ends; ends2];
        
        % sort them
        starts = sort(starts);
        ends = sort(ends);

        
        % how many sections?
        nsec = numel(starts);
        
        slope_pc = nan(nsec,1);
        y1_all = nan(nsec,1);
        y2_all = nan(nsec,1);
        
        
        for i = 1:nsec
            % get the inds
            inds = 1:nvals; inds = inds(:);
            start_ind = inds(x > starts(i)); % gives the index associated with the correct *timestamp*
            end_ind = inds(x < ends(i));
            use_ind = intersect(start_ind,end_ind);
            
            xv = x(use_ind);
            yv = y(use_ind);
            
            % do we have enough points?
            if numel(xv) >= 5
                if strcmp(fit_method,'r')
                    b = robustfit(xv,yv);
                    
                   % find the change along the Y axis relative to intercept
                   y2 = b(2) * xv(end) + b(1); % b(1) is the intercept value
                   y1 = b(2)  *xv(1)   + b(1);
                    
                elseif strcmp(fit_method,'p')
                    b = polyfit(xv,yv,1); % linear 
                    
                    y2 = b(1) * xv(end) + b(2);
                    y1 = b(1) * xv(1)   + b(2);
                end 
                
               y1_all(i) = y1;
               y2_all(i) = y2;
               
               % 2013.09.04 - percentage change from y1 to y2 (the current point)
               slope_pc(i) = 100 * (y2 - y1) / y1; 
            else
                slope_pc(i) = NaN;
            end
       
        end
        
end


% plots

color_it = 1; % use some color
if color_it == 1
    cc = [1 0 0];
elseif color_it == 0
    cc = [0.5 0.5 0.5];
end
if plot_it == 1
    figure(fig_num)
    %subplot(2,1,1)
    if strcmp(split_method,'all')        
        plot(y,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0],'MarkerSize',4,'Marker','o','Color',[0 0 0]) % ignore x-values for simplicity
    else
        plot(x,y,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0],'MarkerSize',4,'Marker','o','Color',[0 0 0]) % we need to preserve TIME
    end
    
    if strcmp(split_method,'all')
        % plot all points on separate plot
        %subplot(2,1,2)
        plot(slope_pc,'Color',cc)
    else
        %subplot(2,1,1)
        hold on
        for i = 1:nsec
            plot([starts(i) ends(i)],[y1_all(i) y2_all(i)],'Color',cc,'LineWidth',1)  
        end
    end
    hold off
end


% get the min and max value
sneg = min(slope_pc);
spos = max(slope_pc);


output.x_start = starts(:);
output.x_end = ends(:);
output.y_start = y1_all(:);
output.y_end = y2_all(:);

output.slope_pc_all = slope_pc; % the individual values
output.slope_pc_neg = sneg;
output.slope_pc_pos = spos;
output.slope_pc_max_abs = max(abs(slope_pc));

if abs(sneg) > abs(spos)
    output.slope_pc_max_sign = sneg;
else
    output.slope_pc_max_sign = spos;
end