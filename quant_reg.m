function out = quant_reg(x,y,keep,xspace,ord,upperP,censor,remap)

% RJE wrapper to perform quantile regression
% makes use of ncquantreg.m from Matlab Central
%
% upperP = vector of percentile values *above the median*, from HIGHEST to LOWEST
%          * note: the default y-axis interpolation size will be 100 - max(upperP)
% censor is a value P; it will cut out x-values below the Pth and (100-P)th percentiles
%
% note: max(upperP) can't be too high if numel(x) is too low and model order is too high; e.g., need at
% least 50000 points to evaluate 99.9th percentile
%
% Specifically, if you see the error "NaN found in Y, interpolation at
% undefined values", you need to either reduce the model order or lower the
% highest upperP values.
%
% note: we *don't* plot in lin or log space; the x-values must themselves
% be in the correct space for the model to work correctly!
%
% RJE, 2016.03.01

x_bins = 100;

% copy of original data
x_orig = x;
y_orig = y;

if isempty(keep)
    keep = ones(numel(x),1);
end

% remove cases we don't want
x(keep == 0) = [];
y(keep == 0) = [];

% also remove NaNs
nanx = isnan(x);
nany = isnan(y);
nanz = sum([nanx nany]) == 2;

x(nanz) = [];
y(nanz) = [];

if isempty(censor)
    % keep all data
else
    minp = prctile(x,censor);
    maxp = prctile(x,100-censor);
    
    keep = x > minp & x < maxp;
    
    x = x(keep);
    y = y(keep);
    
end

if isempty(remap)
    remap = 0;
end

numvals = numel(x);

% 

% check if x-axis is integers
check_int = sum(x == round(x)) / numel(x);

if check_int > .9
    % if 90% of values are integers
    yes_int = 1;
else
    yes_int = 0;
end

% now we can take log of x if we need to
if strcmp(xspace,'log')
    x_orig = log10(x_orig);
    x = log10(x);
end

xmin = min(x);
xmax = max(x);

% now we also need to decide if we want to take log of x
if yes_int == 1
    % just use original values
    x_int  = unique(x);
    x_full = unique(x); 
    
elseif yes_int == 0
    x_int =  linspace(xmin,xmax,x_bins); 
    x_full = linspace(xmin,xmax,1000);
end

num_per_row_temp = numel(x_int);

num_per_row_full = numel(x_full);

%% determine the best order?
if isempty(ord)
   % determine the best low-integer order (1 to 5)
   ord_check = [1 2 3 4 5];
   
   plot_this = 1;
   SS_residuals = polyval_SS(x,y,ord_check,plot_this)
   
  
   ord = input('Select the desired order (1 to 5): ');

end


%% now move on to the desired calculations
% note: ncquantreg.m will break down if there are multiple quantiles chosen
% at the same time and data is too large. So let's iteratively find
% quantiles that are far apart, make sure they don't overlap, and then
% brute force interpolate the rest of them

% manual sequence - just to save time
lowerP = 100 - upperP; % this will be lowest to highest

numP = numel(upperP); 

% figure out the temporary row values
col_temp = [50 upperP lowerP]; 
col_temp = col_temp(:);
col_temp = unique(col_temp); % will sort it as well
num_per_col_temp = numel(col_temp);

% temporary hold
y_curves = nan(num_per_col_temp,num_per_row_temp);

progressbar(0,0)

% get the median
res = ncquantreg_rje(x, y, ord, 50/100); % coefficients are from lowest to highest (i.e., intercept comes first)
res = res(ord+1 : -1 : 1);

med_row = find(50 == col_temp);
y_med = polyval(res,x_int);
y_curves(med_row,:) = y_med;
    
figure(10)
clf
plot(x_orig,y_orig,'.','MarkerSize',4,'Marker','.','LineStyle','none')
hold on
plot(x_int,y_med,'g','LineWidth',1.5)
axis([min(x_orig) max(x_orig) min(y_orig) max(y_orig)])



%% calculate upper and lower target percentile curves
for j = 1:2
    
    if j == 1
        thisP = upperP;
    elseif j == 2
        thisP = lowerP;
    end
    
    for k = 1:numP    
        
        res = ncquantreg_rje(x, y, ord, thisP(k)/100); % coefficients are from lowest to highest (i.e., intercept comes first)
        res = res(ord+1 : -1 : 1);
        
        % evaluate the line
        y_int = polyval(res,x_int);
        
        if k == 1 % starting value
           y_check = y_med; 
        end
        
        % do we keep it?
        y_diff = sign(y_int - y_check);

        if std(y_diff) == 0 % i.e., all the same sign

            % store the result
            this_row = min(find(col_temp >= thisP(k))); % in case of decimal errors
            y_curves(this_row,:) = y_int;

            % redefine the check, and continue
            y_check = y_int;
            
            figure(10)
            plot(x_int,y_int,'r','LineWidth',1)
            axis([min(x_orig) max(x_orig) min(y_orig) max(y_orig)])
        else
            fprintf(['\n   Warning: ' num2str(thisP(k)) 'th percentile overlaps with another percentile.'])
            % keep the same check line and move on to next possible value
        end   

        if j == 1
            progressbar(k/numP,[])
        elseif j == 2
            progressbar([],k/numP)
        end
  
    end % k loop
end % j loop

% remove any rows that are NaN!
nan_row = isnan(y_curves(:,1));

y_curves(nan_row,:) = [];
col_temp(nan_row) = []; % remove this value as well

progressbar(1)

figure(10)
hold off

%% interpolation
% figure out the order
Pmax = max(upperP);
Pmin = 100 - Pmax;
Pdiv = 100 - Pmin; % e.g., 1, 5, 10
col_full = Pmin:Pmin:Pmax;

% also take the calculated values
col_full = [col_full(:); col_temp(:)];
col_full = unique(col_full);

num_per_col_full = numel(col_full);

%%%%%%%% first we do y-axis interpolation
y_curves_temp = nan(num_per_col_full,num_per_row_temp);

% interpolate each column, using calculated values as waypoints
progressbar(0,0)

for j = 1:numel(x_int)
    this_col = y_curves(:,j);

    col_int = interp1(col_temp,this_col,col_full,'linear'); % linear is simplest; also the y-axis space has linear spacing
    
    y_curves_temp(:,j) = col_int;
    progressbar(j/numel(x_int),[])
end

%%%%%%%% now we do x-axis interpolation
y_curves_final = nan(num_per_col_full,num_per_row_full);

for k = 1:num_per_col_full
   this_row = y_curves_temp(k,:);

   y_curves_final(k,:) = interp1(x_int,this_row,x_full,'spline'); % spline is a good choice since we are modeling polynomials
   progressbar([],k/num_per_col_full)
end

progressbar(1)

%% let's actually compute the amount of data below each curve
% note: the code only words when x-values are *discrete* (e.g., integer)
% this will let us know how well the model actually works
% just use the plotted curves for simplicity

if yes_int == 0
    obs_prc_below = 'Not evaluated.';
else

    obs_prc_below = nan(num_per_col_temp,1);
    
    progressbar(0)
    for j = 1:numel(col_temp)
       this_row = col_temp(j) == col_full;
       this_curve = y_curves_final(this_row,:);

       obs_below_temp = nan(1,numel(x_int));

       for k = 1:numel(x_int) % for each unique value

          tar_y = this_curve(k); % the theoretical value 
          obs_y = y(x == x_int(k)); % all observed y-values at this x-value 

          obs_below_temp(k) = sum(obs_y < tar_y);
       end

       obs_below_total = sum(obs_below_temp);

       % turn into percent
       obs_prc_below(j) = 100 * obs_below_total / numel(x);
       progressbar(j/num_per_col_temp)
    end

    progressbar(1)
end


%% remapping the scores
if remap == 1
    progressbar(0)
    fprintf('\n Remapping y-vales ...')

    % set up the final vector
    y_remapped = nan(numel(y_orig),1);

    for i = 1:numel(x_orig)

       % make sure to NaN cases we wish to exclude
       if keep(i) ~= 1
           % final value is NaN
       elseif keep(i) == 1

           this_x = x_orig(i);
           this_y = y_orig(i);

           if this_x < xmin
               % stays at NaN
           elseif this_x > xmax    
               % stays at NaN
           else
              % find the target column
              col_ind = find(x_full >= this_x, 1 ); % find the first (lowest P-value) match

              this_col = y_curves_final(:,col_ind);

              % make sure we are in range
              if this_y > max(this_col)
                  this_prc = 100;
              elseif this_y <= min(this_col)
                  this_prc = min(col_full); % careful: needs to be like this; we can only assign a score based on
                                            % calculated curves; the curve tells us what % of the data
                                            % falls at or below that value
              else 
                  % find the closest percentile value
                  this_prc = min(col_full(this_col >= this_y));         
              end

              if isempty(this_prc)
                  % just a final safety
              else
                  y_remapped(i) = this_prc;
              end
           end
       end 
       progressbar(i/numel(x_orig))
    end

    progressbar(1)
    
    figure(15)
    plot(x_orig,y_remapped,'.','MarkerSize',4,'Marker','.','LineStyle','none')
    
    % correlation
    c = nancorr(x_orig,y_remapped);
    
end % remap option

fprintf('\n\n Finished.\n')

%% output

if isempty(ord)
    out.best_order = ord;
    out.best_R2 = max(R2);
end

out.valid_pts = numel(x);

out.x_plot = x_int';
out.y_plot = y_curves;

out.obs_prc_below = obs_prc_below;

out.x_full = x_full';
out.y_full = y_curves_final;

out.x_orig = x_orig;
out.y_orig = y_orig;

if remap == 1
    out.y_remapped = y_remapped;
    out.Spearman_rho = c.spearman_coef;
end
