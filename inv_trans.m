function [Y, Yrow] = inv_trans(P,N,M,adjust,plot_it)

% [Y, Yrow] = inv_trans(P,N,M,adjust,plot_it)
%
% Use inverse transform sampling to create a matrix of values which follow the *theoretical* PDF specified in P.
% For example, P = [.2 .1 .3 .4] for a 4-level star rating system.
% If there is only one element in P, it will be assumed that this is the binomial sucess rate.
%
% See: https://en.wikipedia.org/wiki/Inverse_transform_sampling
%
% RJE | last edit Nov 2017

if numel(P) == 1
    is_bino = 1;
    P_orig = P;
	P = [1-P_orig P_orig]; % so we have PDF for 0 (failure) and 1 (success)
else
    is_bino = 0;
end

if nargin < 3
    M = 1; % assume we just want an [N x 1] vector
end

if nargin < 4
    adjust = 0; % don't make adjustment (only in binomial case) to make number of 1s exact)
end

if nargin < 5
    plot_it = 0; % don't show plot 
end

% make sure P is normalized from 0 to 1
sumP = sum(P);
P = P ./ sumP;

numP = numel(P);

sumP = cumsum(P);

% assign the values; if only 2 levels; assign 0 and 1

vals = 1:numel(P); % keep this with 1 as minimum for now since it makes life easier later

% fill a matrix of N rows and M columns with numbers from the uniform
X = rand(N,M);

% set up the output variable - fill with NaN, not zeros
Y = nan(N,M);

% do the inverse transform
for p = numP:-1:1 % go backwards
    Y(X <= sumP(p)) = vals(p); % replace values directly in the matrix
end

% now adjust final score
if numP == 2
    % we want 0s and 1s instead
    Y = Y - 1;
else
    % do nothing
end

%% optional correction for simple binomial case
if is_bino == 1 && adjust == 1 && M == 1
    
    % count 1s
    count1 = sum(Y == 1);

    % turn P_orig into a count
    tar_count = round(P_orig * N); % need a whole number since we will replace values directly
        
    % compare to theoretical
    if count1 == tar_count
        % everything is fine
        
    else

        if count1 < tar_count % too few 1s
            
            % turn some 0s into 1s
            change_to = 1;
            
            % get indices of 0s
            ind = find(Y == 0);
            
            % how many need to change?
            nchange = tar_count - count1;
            
        elseif count1 > tar_count % too many 1s

            % turn some 1s into 0s
             change_to = 0;
            
            % get indices of 1s
            ind = find(Y == 1);
            
            % how many need to change?
            nchange = count1 - tar_count;
        end
        
        % do it
        Y(ind(1:nchange)) = change_to;
         
    end
    
    
end

%% optional plot
if plot_it == 1
   figure(10)
   ecdf(Y(:)) % the full matrix
end

%% outputs
Yrow = mean(Y,2); % take the ROW means

