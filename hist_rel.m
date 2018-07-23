function out = hist_rel(data,edges,fignum,draw_box,ymax)

% RJE | 24 Jan 2018

% assume just a single set of data
data = data(:);

minx = min(data);
maxx = max(data);

if nargin < 2
    edges = [];
end

if isempty(edges)
    edges = linspace(minx,maxx,20);
end

if nargin < 3
    fignum = 10;
end

if nargin < 4
    draw_box = 1;
end

if nargin < 5
    ymax = 0; % will make it auto-draw
end

%% calculations

counts = histcounts(data,edges); % histcounts is improved; values that are == final bin are placed *in* the final bin

% normalize it
normc  = counts / numel(data);

% plot it
figure(fignum)
bar(edges(1:end-1),normc);

% change the color map
colormap([.5 .9 .8])

% set figure background to white
set(gcf,'color','w');

xlabel('data')
ylabel('Proportion')

% do this either way
prc = [5 25 50 75 95]; % the percentiles we care about
v   = prctile_nist(data,prc);

% let's add a boxplot on top of this for even more clarity
if draw_box
    hold on
    % figure out where to put it
    if ymax > 0
        c = ymax - .05;
    else
        miny = ceil(max(normc)*20)/20; % divisible by .05
        
        if miny <= .8
            c = miny + .05;
        else
            c = .9;
        end
    end
    
    % width of box
    w = .02;
    
    % adjust yaxis
    ylim([0 c+.05]);
   

    % draw the low whisker
    plot([v(1), v(1)],[c-w/2, c+w/2], 'LineWidth', 1.0,'Color',[0 0 0]);

    % draw the high whisker
    plot([v(5), v(5)],[c-w/2, c+w/2], 'LineWidth', 1.0,'Color',[0 0 0]);

    % draw median
    plot([v(3), v(3)],[c-w, c+w], 'r', 'LineWidth', 1.5,'Color',[1 0 0]); 

    % draw box
    plot([v(2), v(2), v(4), v(4), v(2)], [c-w, c+w, c+w, c-w, c-w], 'LineWidth', 1.0,'LineStyle','-','Color',[0 0 0]);

    % draw lines connecting whiskers to box
    plot([v(1), v(2)],[c, c], 'LineWidth', 1.0,'LineStyle','--','Color',[0 0 0]);
    plot([v(4), v(5)],[c, c], 'LineWidth', 1.0,'LineStyle','--','Color',[0 0 0]);
    
    % highlight first and third quartile
    plot([v(2), v(2)], [c-w, c+w], 'LineWidth', 1.5,'Color',[0 0 1]);
    plot([v(4), v(4)], [c-w, c+w], 'LineWidth', 1.5,'Color',[0 0 1]);
    
    title('Relative histogram + boxplot')
    
    hold off
else
    title('Relative histogram')
end
%% outputs
out.prc    = prc;
out.vals   = v;