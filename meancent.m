function [centered] = meancent(X)

% mean center N x M data (e.g., to validly do a regression)

% can only do an N x 1 or N x M matrix

if numel(size(X)) > 2
   fprintf(' This function only works on [N x 1] or [N x M] data. \n')
   return
end

% we assume we want the mean of the columns (i.e., group average)

mndata = mean(X);

% now just make that many cases
N = size(X,1);
centered = X - (ones(N,1) * mndata);
