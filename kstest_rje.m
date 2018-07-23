function out = kstest_rje(x,dist,plot_it,alpha,type)

% enables other shapes besides normal to be evaluated; we find the
% theoretical values based on observed ECDF values and create a second
% distribution. Advantage: we always get the same KS test result given the
% same x input (i.e., doesn't involve creation of a *random* second sample)

if nargin < 5
    type = 'unequal';
end

if nargin < 4
    alpha = .05;
end

if nargin < 3
    plot_it = 1;
end


x = x(:);
x = sort(x);

% get rid of NaNs
x(isnan(x) == 1) = [];
numx = numel(x);

% get the observed stats
mn_obs = mean(x);
sd_obs = std(x);

% get the *unique* ECDF-values associated with the observed x-values
[ecdf_y x_u] = ecdf(x); % will preserve original x order and frequency (if repeated values)

% we remove 0 and 1, since we can't get theoretical values here
incl = ecdf_y > 0 & ecdf_y < 1; 

x_u(incl == 0) = [];
ecdf_y(incl == 0) = [];

% optional: decimate the data if we think there are too many points


% now find the theoretical x-values at the target y-values

if strcmp(dist,'n')
    x_theo = norminv(ecdf_y);
elseif strcmp(dist,'u')
    x_theo = unifinv(ecdf_y); % this should be the same as ecdf_y, but let's be consistent
else
    % can specify other shapes here
end

% rescale and shift
x_theo = x_theo - mean(x_theo);
x_theo = x_theo / std(x_theo);
x_theo = x_theo * sd_obs; % rescaled
x_theo = x_theo - mean(x_theo); % mean to zero
x_theo = x_theo + mn_obs; % shifted


% now we run a 2-sample KS test
[h pval ks_stat] = kstest2(x,x_theo,alpha,type);

if plot_it == 1
   figure(15)
   plot(x_theo,ecdf_y,'b') % theoretical
   hold on
   plot(x_u,ecdf_y,'r') % observed
   hold off 
end

out.num_x = numx; % useful to count, in case there are NaNs
out.x_mean = mn_obs;
out.x_std  = sd_obs;
out.sig = h;
out.p_val = pval;
out.ks_stat = ks_stat;
out.y_ecdf = ecdf_y(:);
out.x_u = x_u(:);
out.x_theo = x_theo(:);

