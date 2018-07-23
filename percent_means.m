function output = percent_means(data,nan_handling)

% percentage values
p = [.1 .3 .7 .9];

finp = numel(p) + 1;

% how do we deal with nans?

if nargin < 2
    nan_handling = 1;
end

if nan_handling == 0
    % each column will be handled separately
elseif nan_handling == 1
    data_nan = isnan(data);
    sum_nan = sum(data_nan,2); % take the sum across columns
    
    data = data(sum_nan==0,:); % only those rows that don't have NaNs
end

% size of data
nrow = size(data,1);
ncol = size(data,2);

inds(1) = 1;
inds(2:numel(p)+1) = round(p * nrow);
inds(numel(p)+2) = nrow;


output = nan(finp,ncol);

% loop
for c = 1:ncol
   vals = data(:,c);
   vals = sort(vals);
   
   for p = 1:numel(inds)-1
       output(p,c) = mean(vals(inds(p):inds(p+1)));
   end
   
end
    
    