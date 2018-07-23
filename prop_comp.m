%
% A modification of the two-proportion z-test to account for spatially corrolated data.
% Outputs the TWO-tailed p-value associated with the z-test.
% 
% prop_comp(x1,x2,V,k)
%
% x1 = voxels in sample 1
% x2 = voxels in sample 2
%  V = total brain volume in voxels
%  k = expected cluster size
%
% 
% Copyright (C) Rob Ellis | June 2010

% clear variables

%p1 = []; p2 = []; k = []; v = []; x1 = []; x2 = [];

function pval = prop_comp(x1,x2,V,k)

pval=zeros(length(x1),length(x2));

for i=1:length(x1)
    for j=1:length(x2)
      p1 = x1(i)/V;
      p2 = x2(j)/V;

      zval = (p1-p2) / sqrt((p1*(1-p1) + p2*(1-p2)) / (V/k)); 

% now obtain the probability associated with this z-value from the normal
% distribution. multiply value by 2 to get the 2-tailed value, which must
% be less than 0.05 to be significant.

    %pval(i,j)=x1(i)+x1(j); 
    pval(i,j) = 2*(1-normcdf(abs(zval)));
    end
end

 pval;
 zval;

 
end