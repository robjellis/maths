function output = wilcoxon_boot(data,step,iter)

% A novel means to test whether M different repeated measurements of the same
% group of N subjects significantly differ, when different subsets of N are
% taken.
%
% * data is N rows by M columns; NaN cases will be removed across *subjects*
% * group is a vector of numbers (default = 10:10:N)


% get rid of NaNs
data_nan = isnan(data);
data_nan = sum(data_nan,2);

data = data(data_nan==0,:);

% how many subs now?
nsub = size(data,1);
ntest = size(data,2) - 1; % first column is ground truth

group = 10:step:nsub; % auto

% just cap this at 50
group = group(group <= 50);

ngroup = numel(group);

%% loop

% for a given group
matrix = nan(iter,ntest); % iter = rows, tests = columns

% row sum across groups
group_out = nan(ngroup,ntest);
for s = 1:numel(group)  
    for i = 1:iter
        
        % get a random set of subjects
        subs = randperm(nsub);
        subs = subs(1:group(s));
        subs = subs(:);
        
        data_iter = data(subs,:);
        
        for t = 1:ntest
            [p h] = signrank(data_iter(:,1),data_iter(:,t+1));
            matrix(i,t) = h; % will be a 0 or 1 (significant)       
        end
            
    end
    
    % take the SUM of this iter, turn into percent
    group_out(s,:) = 100 * sum(matrix) / iter;
end


%% outputs
output.matrix = matrix;
output.group_out = group_out;
output.group_size = group;
