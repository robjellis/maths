function M = pdf2vals(P,N,M,perm_it,plot_it)

% * Use brute force to create a matrix M of values which follow the *theoretical* PDF specified in P.
%   For example, P = [.2 .1 .3 .4] for a 4-level star rating system with values {1,2,3,4}.
%
% * Setting perm_it = 1 (default) will take a random permutation before reshaping
%   into the final matrix M.
%
% * If there is only one element in P, it will be assumed that this is the
%   binomial sucess rate; resultant values will be either 0 or 1.
%
% RJE | last edit Nov 2017

numP = numel(P);

if nargin == 1 && max(P) > 1
    % we infer N 
    C = P;
    N = sum(C);
end    
    
if nargin < 3
    M = 1; % assume we want a column vector
end

if nargin < 4
    perm_it = 1;
end

if nargin < 5
    plot_it = 0;
end

% total number of values
tot = N*M;

% do we have PDF or counts?
if max(P) < 1
    is_count = 0;
else
    is_count = 1;
end

% special case
if numP == 1 && is_count == 0
    P_orig = P;
	P = [1-P_orig P_orig]; % so we have PDF for 0 (failure) and 1 (success)
    
    % turn PDF into counts
    C = round(P .* tot); % note: may not be exact
    
elseif numP == 1 && is_count == 1
    C = [tot - P, P];
    
elseif numP > 1 && is_count == 0
    % make sure we sum to 1.0
    if sum(P) ~= 1
        error('The sum of probabilities in P does not equal 1.0. Please respecify.')
    end
    
    % turn PDF into counts
    C = round(P .* tot); % note: may not be exact
    
elseif numP > 1 && is_count == 1
    % don't need to do anything except a quick check
    C = P;
    if sum(C) == tot
        % no problems
    else
        error('The sum of counts in P does not match the size [N x M]. Please especify.')
    end
end


%% get the values

% create a large vector 
V = nan(sum(C),1);

ind = 1; % index position

for c = 1:numel(C)
    this_count = C(c);
    this_val = c;
    V(ind:ind+this_count-1) = this_val;
    
    ind = ind + this_count;
end

% special case
if numel(C) == 2
    % assume we want 0s and 1s, not 1s and 2s
    V = V - 1;
end

if perm_it == 1
    % random permutation - optional
    V = V(randperm(numel(V)));
end

% take the first V cases
V = V(1:tot);

% reshape into matrix of [N x M]
M = reshape(V,[N,M]);

%% optional plot
if plot_it == 1
   figure(10)
   ecdf(Vperm) % the full matrix
end