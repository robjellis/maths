function out = smooth_rje(x,y,max_div)

% asummes equal x-axis spacing (i.e., percentage of x-axis space) for simplicity; this allows for lin or
% non-linear "real" spacing
%
% RJE 2016.02.28

% defaults
min_div = 2;
xdiv = 100000;
method = 'linear';



% make sure same size
if numel(x) == numel(y)
    % ok
else
    error('x and y are not the same size')
end
    
numy = numel(y);

% linear interpolation of original y-values
x_orig = linspace(0,100,numy);
x_int  = linspace(0,100,xdiv); % this will make the 4-decomal rounding easier
x_plot = linspace(x(1),x(end),xdiv);

y_int = interp1(x_orig,y,x_int,method);

figure(40)
clf
plot(1:xdiv,y_int,'b')
hold on

% what is the highest number of divisions?
max_div = min(max_div, floor(numy / 2)); % take the smaller value

all_div = min_div:1:max_div;
num_div = numel(all_div);

% index of x-values
case_inds = 1:numy;

% get the full set of possible x-coords (i.e., percentage from start to finish along the x-axis
all_p_cent  = nan(num_div,max_div); % will pad with NaNs but that's fine
all_p_left  = nan(num_div,max_div);
all_p_right = nan(num_div,max_div);

full_y = nan(num_div,xdiv);

progressbar(0)
for j = 1:num_div
    num_bins = all_div(j);
    bin_prop = 1 / num_bins; % e.g., 50 bins = .02 of data
   
    % how many cases per half bin (1M scale)
    half_bin = floor((bin_prop / 2) * xdiv);
    
    % get the bin centers 
    bin_min = bin_prop / 2;
    bin_max = 1 - bin_prop / 2;

    bin_cent = linspace(bin_min,bin_max,num_bins); % these are percentile values  
    
    % scale to 1M 
    case_cent = round(bin_cent * xdiv);
    
    % get the left and right - will already be integers from 0 to 1M
    case_left  = case_cent - half_bin;
    case_right = case_cent + half_bin;
    
    % can't have an index of 0 ...
    case_left(case_left == 0) = 1;
    
    % get the data  
    this_y = nan(num_bins,1);
    
    for k = 1:num_bins
        this_data = y_int(case_left(k):case_right(k));
    
        % save the mean
        this_y(k) = mean(this_data);
    end
    
    figure(40)
    plot(case_cent,this_y,'r')

    % interpolate it
    this_x_int = case_cent(1):case_cent(end);
    this_y_int = interp1(case_cent,this_y,this_x_int,method);
    
    % save it
    full_y(j,this_x_int) = this_y_int;
    
    progressbar(j/num_div)
end

progressbar(1)

% get the mean at each column
full_y_stat = nanmean(full_y,1); 

figure(40)
plot(1:xdiv,full_y_stat,'g')
hold off

return
% round to 4-decimals




for j = 1:num_div
   
    % we don't permute the original data; the order of x and y values is
    % critical, so it stays fixed
    

    
    %num_cent = numel(bin_cent);
    
    % get the bin centers - this is in terms of indexes
    %case_cent_all = prctile(case_inds,bin_cent)


end


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