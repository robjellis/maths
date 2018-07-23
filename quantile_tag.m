function [output] = quantile_tag(data,nquant)

% read in data
% find the desired quantiles (the actual bins)
% tag the data with the bin
% output the bins and tags

%% data check

data = data(:);
ndata = numel(data);

% define quantiles
prc = linspace(0,100,nquant+1);

% get the quantiles

quants = prctile_nist(data,prc);
tags = zeros(ndata,1) + nquant + 1; % start at top

for q = 2:(nquant+1) % nquant+1 is 100 percentile

   tmp = data <= quants(q);
   
   tags = tags - tmp; % will decrement as we march up the quantiles
end


%% outputs

output.quantiles = quants(2:end-1); % ignore the 0th and 100th percentile values
output.tags = tags;