function output = ecdf_manual(data,OPT,plot_it)

% manual ECDF with user-specific params
%
% ** note: only works for "data" that does not have NaNs!
%
% OPT is either (a) [xmin xmax npoints] as a triplet or
%                (b) = 2, which uses regular ecdf.m and then just resorts
%                         the matrix so that values line up with original ordering;
%                         is *much* faster than (d)
%                (c) = 1, which takes unique(data) before doing the ECDF; saves time
%                (d) = 0, which to evaluate ECDF at each value in data; doesn't do unique() 
%                         and thus doesn't sort x-values; this is useful to more easily pair up the original
%                         data values with the resultant ECDF

% just to be safe
data = data(:);


% check for NaNs
if sum(isnan(data))>0
    error('Warning: data contains NaNs and may not sort properly.')
end

% remove NaNs
data(isnan(data)) = []; % still preserve the original data


if nargin < 3
    plot_it = 0;
end

% if nargin < 4 
%     color = [0 0 1];
% else
%     if strcmp(color,'r')
%         color = [1 0 0];
%     elseif strcmp(color,'k')
%         color = [0 0 0];
%     else
%         color = [0 0 1];
%     end
% end

%% which method to do?
if isempty(OPT)
    OPT = 1;
end

if OPT == 0
    x = data;
    npts = numel(x);
elseif OPT == 1;
    x = unique(data);
    npts = numel(x);    
elseif OPT == 2 
    % don't do anything here
elseif numel(OPT) == 3
    npts = OPT(3);
    x = linspace(OPT(1), OPT(2), npts);
    x = x(:);
else
    % no valid method
    return
end

%% loop

if OPT == 0 || OPT == 1
    f = nan(npts,1);

    for i = 1:npts
        f(i) = sum(data <= x(i));
    end

    % divide by size to get scale from 0 to 1
    f = f / numel(data);
    
elseif OPT == 2
    [e x] = ecdf(data);
    
    % x values will be unique(data), but have a double value of the lowest x
    
    e = e(2:end);
    x = x(2:end);
    
    % now we need to match up data with x and sub in e
    % the only way to do this fast is simply by resorting the matrix
    
    % first get the correct indexing
    ind = 1:numel(data);
    ind = ind(:);
    
    % merge it with the values
    mat = [data ind];
    
    % sort it 
    mat = sortrows(mat,1);
    
    % now it should match up exactly with x!
    % confirm this
    if sum(mat(:,1) - x) == 0
        % perfect
    else
        error('Problem with re-sorting of values')
    end
    
    % now add in the e values
    mat = [mat e];
    
    % resort it to original position
    mat = sortrows(mat,2);
    
    % pull out the values
    x = mat(:,1);
    f = mat(:,3);
end


%% figure
if plot_it == 1
    figure(5)
    plot(x,f,'.')
    xlabel('x')
    ylabel('Cumulative probability')
    title(' Empirical CDF')
end

%% output
output.x = x;
output.f = f;
    
    