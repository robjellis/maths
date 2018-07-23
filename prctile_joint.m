function out = prctile_joint(x,y,do_smooth)

% calculate geometric mean and joint percentiles, and show quantile curves
% note: x and y values should have a range 0 to 100
% 
% NaNs in either x or y will result in NaN for geometric mean or joint percentile

if nargin < 3
  do_smooth = 0;
end

curves = (10:10:90)';

% we are either on a 0 to 1 or 0 to 100 scale
if max(x) <= 1
    fprintf('\n Rescaling by a factor of 100 ...')
    x = x * 100;
    y = y * 100;
end

fprintf('\n Calculating geometric mean of x and y...')

G = geomean([x y],2); % if either is NaN, geomean will be NaN
    
    
%% Joint percentile
fprintf('\n Calculating joint percentile ...')

% just do this brute force

nval = numel(x);

nanx = isnan(x);
nany = isnan(y);

nanz = sum([nanx nany],2) == 0; % = 1 only if both are not NaN

J = nan(nval,1);

progressbar(0)
for i = 1:nval
    progressbar(i/nval)
    if nanz(i) == 1 % means good to go

        thisx = x(i);
        thisy = y(i);

        J(i) = 100 * sum(sum([x<thisx y<thisy],2)== 2) / nval;
    else 
        % ... exclude the case
    end
    
end

progressbar(1)

%% calculate curves
% let's round to make life easier

if nval < 50000
    % round to whole number
    mult = 1;
else
    % round to nearest 0.1
    mult = 10;
end

rx = round(x*mult)/mult;
ry = round(y*mult)/mult;
rG = round(G*mult)/mult;
rJ = round(J*mult)/mult;

figure(301)
clf
plot(x,y,'MarkerSize',4,'Marker','.','LineStyle','none', 'Color',[.3 .3 .3]) % for G curves
hold on

figure(302)
clf
plot(x,y,'MarkerSize',4,'Marker','.','LineStyle','none', 'Color',[.3 .3 .3]) % for J curves
hold on

for k = 1:numel(curves)
    P = curves(k);
    
    % get the data for this curve
    Gx = rx(rG == P);
    Gy = ry(rG == P); 
    
    Jx = rx(rJ == P);
    Jy = ry(rJ == P);
    
    % sort
    Gvec = [Gx Gy];
    Jvec = [Jx Jy];

    Gvec = sortrows(Gvec,1);
    Jvec = sortrows(Jvec,1);
    
 %   do_smooth = 1;
    
    if do_smooth == 1 % simple smoothing on x and y; results look better with more original data
        Gx_sm = smooth(Gvec(:,1));
        Gy_sm = smooth(Gvec(:,2));
        
        Jx_sm = smooth(Jvec(:,1));
        Jy_sm = smooth(Jvec(:,2));
    end
        
    figure(301)
    if do_smooth == 1
        plot(Gx_sm,Gy_sm,'k','LineWidth',1.1)
    end
    plot(Gx,Gy,'k.','MarkerSize',7)
    
    
    figure(302)    
    if do_smooth == 1
        plot(Jx_sm,Jy_sm,'k','LineWidth',1.1)
    end
    plot(Jx,Jy,'MarkerSize',10,'Marker','*','LineStyle','none', 'Color',[1 0 0])
    
end

figure(301)
title('GM deciles')
xlabel('x-values')
ylabel('y-values')
axis([0 100 0 100])
hold off

figure(302)
title('Joint deciles')
xlabel('x-values')
ylabel('y-values')
axis([0 100 0 100])
hold off

figure(303)
plot(G,J,'MarkerSize',4,'Marker','.','LineStyle','none')
hold on
xlabel('GM')
ylabel('Joint percentile')
axis([0 100 0 100])

% plot the theoretical line
Gtheo = 0:100;
Jtheo = Gtheo.^2 / 100;
plot(Gtheo,Jtheo,'r','LineWidth',1.5)
hold off

figure(304)
ecdf(G)
hold on
ecdf(J)
hold off
legend('show','Location','SouthEast')


out.x = x;
out.y = y;
out.G = G;
out.J = J;