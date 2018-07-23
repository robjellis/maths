function [r] = qqplot_rje(x1,x2,nquant,fignum,data1_name,data2_name)

% a function that doesn't require matlab stats toolbox
%
% note: percentiles range from 0 to 100, so fewer decimals are required
% 
% if "fignum" = 0, then no new figure will be called (e.g., can be called
% as a subplot)

stats_ok = which('normcdf');

if sum(size(stats_ok)) > 0 % only if we have the stats toolbox ...

   % generate a series of equal step percentile values based on z scores
   
    z = linspace(-3.5,3.5,nquant);
    nz = normcdf(z);                  % to get corresponding CDF values from 0.0 to 1.0
    p = round(nz*100000);              % round 5 decimal places
    p = unique(p);                    % get rid of repeats
    p = p / 1000;                      % to get back to percentiles

        if p(1) == 0
           p(1) = [];
        end

        if p(end) == 100;
           p(end) = [];
        end

   percentiles = p;
   
else
   % 50 percentile values evenly spaced along Z 
   percentiles = [0.0230000000000000,0.0390000000000000,0.0650000000000000,0.107000000000000,0.170000000000000,0.267000000000000,0.411000000000000,0.621000000000000,0.921000000000000,1.34000000000000,1.91600000000000,2.68900000000000,3.70700000000000,5.02100000000000,6.68100000000000,8.73700000000000,11.2320000000000,14.1990000000000,17.6560000000000,21.6020000000000,26.0160000000000,30.8540000000000,36.0490000000000,41.5160000000000,47.1530000000000,52.8470000000000,58.4840000000000,63.9510000000000,69.1460000000000,73.9840000000000,78.3980000000000,82.3440000000000,85.8010000000000,88.7680000000000,91.2630000000000,93.3190000000000,94.9790000000000,96.2930000000000,97.3110000000000,98.0840000000000,98.6600000000000,99.0790000000000,99.3790000000000,99.5890000000000,99.7330000000000,99.8300000000000,99.8930000000000,99.9350000000000,99.9610000000000,99.9770000000000];

   
end

q1 = prctile_nist(x1,percentiles);
q2 = prctile_nist(x2,percentiles);

if fignum > 0
    figure(fignum)
else
    %
end

plot(q1,q2,'*')
hold on

% now plot the reference line (code borrowed from Matlab qqplot)

q1x = prctile_nist(x1,25);
q3x = prctile_nist(x1,75);
q1y = prctile_nist(x2,25);
q3y = prctile_nist(x2,75);
qx = [q1x; q3x];
qy = [q1y; q3y];

dx = q3x - q1x;
dy = q3y - q1y;
slope = dy./dx;
centerx = (q1x + q3x)/2;
centery = (q1y + q3y)/2;
maxx = max(q1);
minx = min(q1);
maxy = centery + slope.*(maxx - centerx); % passes through center
miny = centery - slope.*(centerx - minx);

mx = [minx; maxx];
my = [miny; maxy];

plot(mx,my,'Color','r','LineStyle','-','Marker','none');
%plot(qx,qy,'Color','r','LineStyle','-','Marker','none');
hold off

% what is the Pearson correlation between the quantiles?
% if we have "normcdf", then we should have "corr"

if sum(size(stats_ok)) > 0
   r = corr(q1,q2);
else 
   r = NaN; 
end


lx = xlabel([data1_name ' quantiles']); ly = ylabel([data2_name ' quantiles']); 
titl = title(['Quantile-Quantile plot' sprintf('\n(') num2str(nquant) ' quantiles)']);
set(lx,'Interpreter','none'); set(ly,'Interpreter','none'); set(titl,'Interpreter','none'); 


