function output = fwer(C,cov,n)

if nargin < 3
    n = 5000000; % rje: this needs to be large (>= 5M)
end

% cov matrix
cmat = eye(C,C);

% replace 0s with cov
cmat(cmat==0) = cov;

% create the multivariate normal 
means = zeros(1,C);

vals = mvnrnd(means,cmat,n);

% we want *exactly* 5% to be significant
prc = prctile_nist(vals,95);

% now threshold each column
sig = zeros(n,C);

for c = 1:C
    sig(:,c) = vals(:,c) > prc(c);
end

fpr = sum(sig) / n; % each value should be exactly .0500

% now we evaluate the FWER for 1:C

fwer_obs = nan(1,C);

for c = 1:C
   fwer_obs(c) = sum(max(sig(:,1:c),[],2)) / n; 
end

% what is the theoretical worst case scenario inflation?
cval = 1:C;
fwer_theo = 1 - (1 - .05).^cval;

output.fwer_obs = fwer_obs;
output.fwer_theo = fwer_theo;