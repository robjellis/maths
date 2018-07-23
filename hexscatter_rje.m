function params = hexscatter_rje(X,Y,do_log2,res,fignum,xlim,ylim)

% RJE wrapper for hexscatter.m (https://github.com/brazilbean/bean-matlab-toolkit)
% to make slightly better plots
%
% calls a modified version of hexscatter.m: hexscatter_mod.m
%
% RJE | 6 Feb 2018

X = X(:);
Y = Y(:);

numX = numel(X);
numY = numel(Y);

if numX == numY
	% ok
else
	error('X and Y vectors are not the same length.')
end

if nargin < 3
	do_log2 = 'true';
end

if nargin < 4
	% set the grid such that res^2 = number of points
	res = ceil(sqrt(numX));
end

if res > 200
	res = 200;
end

if nargin < 5
	figure(100)
else
	if fignum == 0
		% we call the figure from outside the function
	else
		figure(fignum)
	end
end

if nargin < 6
	xlim = [min(X) max(X)];
end

if nargin < 7
	ylim = [min(Y) max(Y)];
end

%% plot
% critical switch!
if strcmp(do_log2,'true')	
	[h, params] = hexscatter_log2(X,Y,'ShowZeros','true','res',res,'xlim',xlim,'ylim',ylim);
else
	% use the unaltered code
	h = hexscatter(X,Y,'ShowZeros','true','res',res,'xlim',xlim,'ylim',ylim);
end

% make it pretty
set(gcf,'color','w');
xlabel('X values')
ylabel('Y values')
colorbar; % show the colorbar

if strcmp(do_log2,'true')
	title('Heat map: log2(counts + 1)')
else
	title('Heat map: raw counts')
end



