function [y x] = histc_gen(X,log_it,plot_type,num_max,num_min,height_min,fig_num,cur_plot,num_plot)

% modification of histc with more options for histograms
% 
% histc_gen(X,log_it,plot_type,num_max,num_min,height_min,fig_num,cur_plot,num_plot)
%
% "plot_type": 'b' for bar, 'l' for line, or 'n' for none
%
% "num_max": max number of bins; should be an EVEN number
%
% "num_min" should be an EVEN number, which will center a bin around median
% (assuming normal distribution)
%
% "bin_min" is a percentage value: how many values (as a percent of total)
% should the minimum bin have? (default = 1.0)
%
% "log_it": 0 for don't, 1 for yes; will plot log-spaced values on y-axis
% but using the 'real' numbers; e.g., 0 1 2 4 8 16 ...
%
% "fig_num": 0 means don't make a new figure; -1 means don't plot at all
%
%
% Rob Ellis | 2013.05.06

X = X(:);

% get rid of NaN
X(isnan(X)) = [];
nX = numel(X);

% we make bins that are sensitive to integers if we determine that the data
% has mostly integer values

numint = sum(mod(X,1) == 0);

if numint / nX > .50
    int_val = 1;
else
    int_val = 0;
end

if nargin < 2
   log_it = 1;
end

if nargin < 3
    plot_type = 'b'; % 'b' = bar, 'l' = line
end

if nargin < 4
   num_max = 50; 
end

if nargin < 5
   num_min = 6;
end

if nargin < 6
   height_min = 1.0; % at least 1.0 percent of total data in least frequent bin
end

if nargin < 7
   fig_num = 0; % will call Figure 1 if no figure is active
end

if nargin < 8
   cur_plot = 1;
end

if nargin < 9
   num_plot = 1;
end

% set colors for hist
cc = hsv(num_plot);


%% we automatically calculate the number of bins
     
% RJE solution to have balanced bins (around the median)
%val = max(max(X)-median(X),median(X)-min(X));
%minx = median(X) - val;
%maxx = median(X) + val; % before, added 1e-10 to this, but this will cause problems if we have whole number integers and are trying to count them!

minx = floor(min(X));
maxx = ceil(max(X));


if num_max == num_min
    % force a single run with this many bins
       
    if int_val == 1
        
        % OK, we have to be more specific about doing bins here
        nint = numel(minx:1:maxx); %
        
        if nint == 1
            fprintf('\n Warning: A variable has no variance. \n\n')
            return
        end
        m = (num_max - nint) / (nint-1); % how many more divisions *between* existing integer values?
        
        mult = [    1       3      9       19      49      99   199];
        val  = [    .5      .25    .1      .05     .02     .01  .005];
        m = ceil(m); % round up
        
        mind = min(find(mult >= m));

        x = minx:val(mind):maxx; % will ensure that we use a division that preserves integers
        
        % shift the edges of each bin half a bin down so that we can center around integers
        xdiff = (x(2) - x(1))/2;
        
        xadj = x - xdiff; 
        xadj(end+1) = x(end) + xdiff; % so we still capture the original value!
        
        x = xadj;
        
    elseif int_val == 0
        % much easier
        x = linspace(minx,maxx,num_max);
        %x(end) = x(end) + 1;
    end
    
        % now do the histc
        y = histc(X,x);

        % get rid of all bins where y = 0
        %x = x(y > 0);
        %y = y(y > 0);

    
else
    for j = num_max:-2:num_min % even number of bins centers a bin around median, assuming normal dist.

        x = linspace(minx,maxx,j); % starting condition only
        x(end) = x(end) + 1; % so that count for x(end) is always 0
        y = histc(X,x);

        % get rid of all bins where y = 0
        %x = x(y > 0);
        %y = y(y > 0);

        if (min(y)/nX)*100 < height_min;
            % use fewer bins
        elseif min(y) == 1 
            % if 1 is the lowest count, then we have too many bins
        else
            % use the *previous* bin division
            break
        end

    end
end

%% now we do the optional log2 transformation
if log_it == 1
    % we do a binary (base 2) log transform, and then some special tricks
    % so that it displays the original units on the y-axis
    
    log_disp = [0    1   2	4	8	16	32	64	128	256     512     1024	2048	4096	8192	16384	32768	65536	131072	262144	524288	1048576  2097152   ];
    log_val  = [-1   0   1	2	3	4	5	6	7	8       9       10      11      12      13      14      15      16      17      18      19      20       21   ];
  
    y = log2(y); % log2(0) will = -Inf and will thus not show up on the plot!
    
 
end

% now plot the result using a bar graph
if fig_num >= 0
    
    if fig_num == 0
        % current figure
    else
        figure(fig_num)
    end
   

    % have to shift the bins *up* by half a bar to display as edges correctly
    % need to do this for both line and bar plots, even with integer data
    xdiff = (x(2) - x(1)) / 2;
    xplot = x + xdiff;


    % what is the highest and lowest counts, and round?
    miny = floor(min(y)); % we end up not using this when doing a log plot, since we force the min log value to display as 0
    maxy = ceil(max(y));

        % make the plot
        if strcmp(plot_type,'b')
            if int_val == 1
                width = 0.8; % then we want arbitrarily small bars!
            else
                width = 1.0;
            end 
                if log_it == 0
                    bar(xplot,y,width,'BaseValue',0);
                elseif log_it == 1
                    bar(xplot,y,width,'BaseValue',-1); 
                end
                
        elseif strcmp(plot_type,'l')
            plot(xplot,y,'LineWidth',1,'color',cc(cur_plot,:),'DisplayName',['Data ' num2str(cur_plot)]);
        else
            % don't plot
        end

    if log_it == 1
        
        % since we are on a log plot, makes sense to always have 0 as the
        % y-minimum value, and also easier to display values = 1
        
        % set the figure tick *limits* like this
        if min(X) == 0
            minxplot = 0;
        elseif min(X) > 0 && min(X) <= 1
            minxplot = 0; % useful for percents with lower value = 0
        else
            minxplot = floor(min(xplot));
        end
        
        if max(X) <= 1
           maxxplot = 1;
        elseif max(X) > 1 && max(X) < 100
           maxxplot = ceil(max(X));
        elseif max(X) <= 100 % useful for percents from 0 to 100
           maxxplot = 100;
        else
           maxxplot = ceil(max(xplot));
        end
        
        axis([minxplot,maxxplot,-1,maxy])
        
        % now set the tick *spacing*, based upon how many ticks we have
        nticks = numel(-1:maxy);
        
        if nticks <=6 % should display fine on a small plot
           yticks = -1:1:maxy;
        elseif mod(nticks,2) == 0 % even
           yticks = -1:2:(maxy+1);
        elseif mod(nticks,2) == 1 % odd
           yticks = -1:2:maxy;
        end

        set(gca,'YTick',yticks) % actual y values

        % redo the labels        
        set(gca,'YTickLabel',log_disp(yticks+2)); % + 2 will do the proper conversion
    elseif log_it == 0
        % specify the limits
        
        %ylab = linspace(miny,maxy,6);
        %set(gca,'YTick',ylab);
        %set(gca,'YTickLabel',ylab);
    end
    
        ylabel('Count')

end

% just to be sure
y = y(:);
x = x(:);