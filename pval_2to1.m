function [out] = pval_2to1(stat,pval,direc)

% convert 2-tailed p-values to 1-tailed p-values, based on the
% expected direction of the result.
% direct: 1 for positive and -1 for negative
% make sure data is correct

if size(stat,1) == size(pval,1)
    % OK
else
    return
end

%

ntests = size(stat,1);
out = zeros(ntests,1);

for i = 1:ntests
    if sign(stat(i)) == direc
       % this is a test we expect
       out(i) = pval(i)/2;
       
    else
       % this is a test we don't expect 
       out(i) = 1 - pval(i)/2; 
    end
    
end