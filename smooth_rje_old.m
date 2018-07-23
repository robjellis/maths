function out = smooth_rje_old(x,y,npts,iter)

% figure out x-range
xmin = min(x);
xmax = max(x);

numx = numel(x);
xint = linspace(xmin,xmax,100);

yint_all = nan(iter,numel(xint));
for i = 1:iter
    
    % get the points
    pts = randint(npts,1,[2 numx-1]);
    
    pts2 = [1 pts' numx];
    pts2 = unique(pts2); % this will be sorted
    
    x_this = x(pts2);
    y_this = y(pts2);
    
    % iterpolate
    yint_all(i,:) = interp1(x_this,y_this,xint,'linear'); 
   
end

% take the average at each point
yint_final = mean(yint_all);

out.y_smooth = yint_final;

figure(300)
plot(x,y)
hold on
plot(xint,yint_final,'r')
hold off