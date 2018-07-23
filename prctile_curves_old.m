function out = prctile_curves_old(x,y,keep_index,bin_space,num_bins,bin_size,boot_samp,boot_iter,smooth_meth,smooth_param)

% prctile_boot(x,y,keep_index,bin_space,num_bins,boot_samp,boot_iter,smooth_meth,smooth_param)
%
% bin_space is either 'case' - equal division among cases, ignore x-values
%                     'lin' 
%                     'log' 
% bin_size is [] if 'case'; otherwise, it is the percentage of total data in each bin
%
% num_bins = how many bins across the x-axis space (default = 100)
% boot_samp = 1000000 by default
% smooth_meth = 'ksr' - kernel smoothing; 
%               'rlowess' - robust linear;
%               'rloess' - robust quadratic; 
%               'moving' - moving average
% smooth_param = for rloess or rlowess, a number from 0 to 1; proportion of data to use at each point
%                for ksr, can leave at [] for auto
%                for 'moving', the number of points (can make it proportional to num_bins)

%% defaults

prc = (.1 : .1 : 99.9)';
prc_plot = [.1 1 10 25 50 75 90 99 99.9]'; 

num_prc = numel(prc);

%% checks
if numel(x) ~= numel(y)
    error('Mismatched number of values')
else 
    % things are fine
end

% calculate bin size
if bin_size * 2 <= num_bins
    % ok
else
    error('Bin size must be <= half the number of bins!')
end
    
%% get rid of cases we don't want

% make a copy, so we can preserve the original order
x_orig = x;
y_orig = y;

% a priori exclusion
if isempty(keep_index)
    keep_index = ones(numel(x),1);
end

x(keep_index == 0) = [];
y(keep_index == 0) = [];

% also get rid of NaNs in x or y
nan1 = isnan(x);
nan2 = isnan(y);
nan3 = sum([nan1 nan2],2) >= 1;

x = x(nan3 == 0);
y = y(nan3 == 0);

min_obs = min(x);
max_obs = max(x);

numx = numel(x);

% check if x-axis is integers
check_int = sum(x == round(x)) / numel(x);

if check_int > .9
    % if 90% of values are integers
    yes_int = 1;
else
    yes_int = 0;
end

% important: permute the vector and then resort by x-values
randvec = randperm(numx);
x = x(randvec);
y = y(randvec);

% now sort by x
zvals = [x y];
zvals = sortrows(zvals,1);

x = zvals(:,1); % now we can easily pull cases
y = zvals(:,2);


% locations of bin centers, in terms of cases
case_inds = 1:numx;

% we either define the locations we look at by equal spacing along cases,
% or by equal spacing along x-values directly; for now, we just do the
% former since it is easier to implement (and to describe)

if strcmp(bin_space,'case')
    bin_prop = 1 / num_bins; % e.g., 50 bins = .02 of data
    half_cases = floor((numx * bin_prop) / 2) - 1;
else
    half_cases = floor((numx * bin_size/100)/2) - 1;
end

% find the first and last center
case_cent_first = half_cases + 1; % since we can't start with index 0
case_cent_last  = numx - half_cases;

% get the x-values at these indexes
xmin = x(case_cent_first);
xmax = x(case_cent_last);
    
if strcmp(bin_space,'case')
    case_cent_all = floor(linspace(case_cent_first,case_cent_last,num_bins));
    x_cent_all    = x(case_cent_all); % these will be used for plotting

else   
    if strcmp(bin_space,'lin')
        x_cent_all = linspace(xmin,xmax,num_bins);      
    elseif strcmp(bin_space,'log')
        x_cent_all = logspace_rje(xmin,xmax,num_bins);
    end 
    
    if yes_int == 1 
        x_cent_all = round(x_cent_all);
        
        % careful! we have to have *distinct* x-axis values for later interpolation steps
        x_cent_all = unique(x_cent_all);
        
        % return a warning if needed
        if numel(x_cent_all) < num_bins
            msg = [' Only ' num2str(numel(x_cent_all)) ' bins are possible between xmin and xmax.'];
            warning(msg)
        end
        
        num_bins = numel(x_cent_all);
        
        
    end
   
    % now need to find the corresponding *case* locations
    case_cent_all = nan(num_bins,1);
    
    for j = 2:num_bins-1
        
       % find the midpoint of cases that have this n; note they are already sorted from lowest to highest

       % careful!: need to find the closest *actual* x-value in order to match
       x_u = unique(x);
       d = abs(x_u - x_cent_all(j));
       find_this_x = x_u(d == min(d));

       case_cent_all(j) = floor(median(find(x == find_this_x)));
       
    end
    
    % manually do the first and last case; otherwise, we may end up with
    % negative cases (or numbers outside actual cases)
    case_cent_all(1) = case_cent_first;
    case_cent_all(end) = case_cent_last;
end

case_min_all = case_cent_all - half_cases;
case_max_all = case_cent_all + half_cases;

% to store the x-value min and max for each bin, just so we can confirm
% whether bins overlap (if they should, depending on the method we use)
x_min_all  = nan(num_bins,1);
x_max_all  = nan(num_bins,1);

items = nan(num_bins,1);
vals  = nan(num_prc,num_bins);

%% 1. do the bootstrap
fprintf('\n Performing the bootstrap percentile calculation ...')
progressbar(0,0,0,0,0)

for i = 1:num_bins
       
    % get the y-values
    this_case_min = case_min_all(i);
    this_case_max = case_max_all(i);
    
    y_set = y(this_case_min:this_case_max);
    
    items(i) = numel(y_set);
       
    % get the x-range
    x_min_all(i) = x(this_case_min);
    x_max_all(i) = x(this_case_max);
    
      
    
    %% now the CI calculation
    
    % this needs to iterate so we can get stable percentile values,
    % especially in the taisl
    
    if isempty(boot_iter)
        boot_iter = 1; % just so the code works
    end
    
    vals_temp = nan(num_prc,boot_iter); % reset for each bin
    
    for k = 1:boot_iter
        
        boot_ind = randint(boot_samp,1,[1 items(i)]);

        % get the distrubtion
        stat_boot = y_set(boot_ind);

        % get the percentiles
        vals_temp(:,k) = prctile(stat_boot,prc);
        progressbar([],k/boot_iter,[],[],[])
    end
    
    % now we take the median of each *row*
    vals(:,i) = median(vals_temp,2); 
    
    progressbar(i/num_bins,[],[],[],[])
end


%% 2. calculate smooth
fprintf('\n Smoothing values ...')

sm_vals      = nan(num_prc,num_bins);

sm_param   = nan(num_prc,1);

for j = 1:num_prc % each percentile is done separately
    this_vals = vals(j,:);
    
    if strcmp(smooth_meth,'ksr')
        h = smooth_param; % can be left as []
        n = numel(x_cent_all);
        
%         if strcmp(bin_space,'lin')
%             sm_x = x_cent_all; % use the original values
%         elseif strcmp(bin_space,'log')
%             sm_x = (1:numel(x_cent_all))'; % works better this way; we will still plot with x_cent_all however
%         end
%         
%         ksr_res = ksr_rje(sm_x,this_vals,h,n);
        
        % rje: ksr works best if x is just ordinal
        ksr_res = ksr_rje(1:numel(this_vals),this_vals,h,n);
        
        sm_vals(j,:) = ksr_res.f;
        sm_param(j)  = ksr_res.h;
    elseif strcmp(smooth_meth,'equal') % all weights equal
        weight_val = 1 / smooth_param;
        weights = ones(1,smooth_param) * weight_val;
        
        sm_vals(j,:) = conv(this_vals,weights,'same');
    
    elseif strcmp(smooth_meth,'quad')
        % need to build this filter
    elseif strcmp(smooth_meth,'moving')
        % let's make some adjustments to the sequence; pad the ends so they
        % are less biased when we smooth it
        pad = floor(smooth_param / 2); % assuming the value is odd
        this_vals_adj = nan(pad + numel(this_vals) + pad, 1);
        
        this_vals_adj(1:pad) = this_vals(1);
        this_vals_adj(pad+1:end-pad) = this_vals;
        this_vals_adj(end-pad+1:end) = this_vals(end);
        
        vals_out = smooth(this_vals_adj,smooth_param,smooth_meth);
        vals_out_trim = vals_out(pad+1 : end-pad);
        sm_vals(j,:) = vals_out_trim;
        
    else
        % Matlab functions:  SMOOTH(X,Y,SPAN,METHOD)
        % rje: for simplicity, ignore x-values; just use ordinal
        sm_vals(j,:) = smooth(1:numel(this_vals),this_vals,smooth_param,smooth_meth); 
    end
    
    progressbar([],[],j/num_prc,[],[])
end

% get the ones to plot

ind = nan(numel(prc_plot),1);

for k = 1:numel(prc_plot)
   ind(k) = find(prc >= prc_plot(k),1); % do this just in case we don't match decimals exactly  
end

vals_plot = vals(ind,:);
sm_vals_plot = sm_vals(ind,:);

figure(200)
if strcmp(bin_space,'lin')
    plot(x_cent_all,x_min_all,'.')
    hold on
    plot(x_cent_all,x_max_all,'.')
    hold off
elseif strcmp(bin_space,'log')
    loglog(x_cent_all,x_min_all,'.')
    hold on
    loglog(x_cent_all,x_max_all,'.')
    hold off
end


% show key prctiles
figure(201)
if strcmp(bin_space,'lin')
    plot(x,y,'MarkerSize',4,'Marker','.','LineStyle','none')
elseif strcmp(bin_space,'log')
    semilogx(x,y,'MarkerSize',4,'Marker','.','LineStyle','none')
end
    
hold on
for k = 1:numel(prc_plot)
    if strcmp(bin_space,'lin')
        plot(x_cent_all,vals_plot(k,:),'m','LineWidth',1)
    elseif strcmp(bin_space,'log')
        semilogx(x_cent_all,vals_plot(k,:),'m','LineWidth',1)
    end
end
hold off


% now plot with rloess results
figure(202)
if strcmp(bin_space,'lin')
    plot(x,y,'MarkerSize',4,'Marker','.','LineStyle','none')
elseif strcmp(bin_space,'log')
    semilogx(x,y,'MarkerSize',4,'Marker','.','LineStyle','none')
end

hold on
for k = 1:numel(prc_plot)
    if strcmp(bin_space,'lin')
        plot(x_cent_all,sm_vals_plot(k,:),'r','LineWidth',1)
    elseif strcmp(bin_space,'log')
        semilogx(x_cent_all,sm_vals_plot(k,:),'r','LineWidth',1)
    end  
end
hold off


%% 3. interpolation
fprintf('\n Interpolating percentile curves ...')

if yes_int == 0
    % just use these x-values
    x_unique = unique(x); % cut out the duplicates

    % only look at "valid" range (xmin and xmax exclude the tails)
    x_unique(x_unique < xmin) = [];
    x_unique(x_unique > xmax) = [];
    
    x_interp = x_unique;
else
    % add in missing inteters
    x_interp = xmin:xmax;
end

y_interp = nan(num_prc,numel(x_interp));

for j = 1:num_prc

   y_interp(j,:) = interp1(x_cent_all,sm_vals(j,:),x_interp,'linear'); % note: we could select 'cubic' to give us pchip 

   progressbar([],[],[],j/num_prc,[])
end


%% 4. remapping the scores
fprintf('\n Remapping y-values ...')

% set up the final vector
y_remapped = nan(numel(y_orig),1);

for i = 1:numel(x_orig)
    
   % make sure to NaN cases we wish to exclude
   if keep_index(i) == 0
       % final value is NaN
   else
       
       this_x = x_orig(i);
       this_y = y_orig(i);

       if this_x < xmin
           % stays at NaN
       elseif this_x > xmax    
           % stays at NaN
       else
          % find the target column
          col_ind = x_interp == this_x; % there will be only one match

          this_col = y_interp(:,col_ind);

          % make sure we are in range
          if this_y > max(this_col)
              this_prc = 100;
          elseif this_y <= min(this_col)
              this_prc = min(prc);  % simple rule is assign the percentile curve that is >= value
                                    % we don't assign a value of zero,
                                    % since nothing can be below the 0th percentile
          else
              this_prc = min(prc(this_col >= this_y));         
          end

          if isempty(this_prc)
              % just a final safety
          else
              y_remapped(i) = this_prc;
          end
       end
   end 
   progressbar([],[],[],[],i/numel(x))
end

progressbar(1)

%% 5. calculate the "simple" percenilte directly from ecdf_mod
y_ecdf_prc = ecdf_mod(y_orig) * 100;

figure(203)
plot(y_ecdf_prc,y_remapped,'MarkerSize',4,'Marker','.','LineStyle','none')
xlabel('Percentile from ECDF')
ylabel('Percentile from bootstrap')

figure(204)
if strcmp(bin_space,'lin')
    plot(x_orig,y_remapped,'MarkerSize',4,'Marker','.','LineStyle','none')
elseif strcmp(bin_space,'log')
    semilogx(x_orig,y_remapped,'MarkerSize',4,'Marker','.','LineStyle','none')
end 

%% 6. compute the amount of data that falls below the target curves
% this is easy now that we have remapped the scores!

nplot = numel(prc_plot);

obs_below_plot = nan(nplot,1);

num_non_nan = sum(isnan(y_remapped) == 0); % we only count the *non* NaN values

for k = 1:nplot
    obs_below_plot(k) = 100* sum(y_remapped <= prc_plot(k)) / num_non_nan;
end

%% correlations
corr1 = nancorr(x,y); % only the data we care about; not the *full* data
corr2 = nancorr(x,y_remapped);
corr3 = nancorr(y_ecdf_prc,y_remapped);

%% outputs
out.obs_min     = min_obs;
out.obs_max     = max_obs;

out.use_min     = xmin;
out.use_max     = xmax;

out.x_cent_all  = x_cent_all;
out.case_cent_all = case_cent_all;

out.sm_param    = sm_param;

out.items   = items;
out.vals    = vals;
out.sm_vals = sm_vals;

out.x_interp = x_interp;
out.y_interp = y_interp;

out.x_orig = x_orig; % just a copy to be safe
out.y_orig = y_orig;
out.y_remapped = y_remapped; % the final bootstrap-adjusted percentile score

out.obs_below_plot = [prc_plot(:) obs_below_plot];
out.spearman_x_y = corr1.spearman_coef;
out.spearman_x_y_remapped = corr2.spearman_coef;
out.spearman_ecdf_y_remapped = corr3.spearman_coef;
out.num_NaN_prc = sum(isnan(y_remapped));


