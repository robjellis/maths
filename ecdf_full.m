function output = ecdf_full(data,min,max,nstep,space,plot_it)

if nargin < 3
    % assume normal data
    min = -4;
    max = 4;
    nstep = 1000; % RJE: if this gets too small, we can't do the simultaneous KS tests later
    space = 'lin';
end

if nargin < 6
    plot_it = 0;
end

% just to be safe
data = data(:);


% check for NaNs
if sum(isnan(data))>0
    error('Warning: data contains NaNs and may not sort properly.')
end

% remove NaNs
data(isnan(data)) = []; % still preserve the original data

x = unique(data);

nx = numel(x);

f = nan(nx,1);

for i = 1:nx
    f(i) = sum(data <= x(i));
end
    

% divide by size to get scale from 0 to 1
f = f / numel(data);


% now we adjust the values

if strcmp(space,'lin')
    xfull = linspace(min,max,nstep);
elseif strcmp(space,'geo')
    xfull = geospace(min,max,nstep);
end

xfull = xfull(:);

nxfull = numel(xfull);

ffull = zeros(nxfull,1);

for i = 1:nx
    
    % indexes we care about
    inds = xfull > x(i);
    
    ffull(inds) = f(i); % this will iterate correctly
    
end

output.x = x;
output.f = f;
output.x_full = xfull; % a column
output.f_full = ffull; % a column

if plot_it == 1
    figure(10)
    plot(xfull,ffull,'k')
    hold on
    plot(x,f,'r.')
    hold off
end
