function [wt p dfe dfer] = welch_t(b1,b2,se1,se2,df1,df2)

% ** note: this is designed for tests with stats coming out of a (robust) 
%    regression analysis, not a Welch test on vectors (MATLAB implements
%    the latter already)
%
% b  = regression intercept (mean)
% se = regression intercept standard error

% need to get a matrix for df values
ncol = size(b1,2);
ddf1 = ones(size(b1));
ddf2 = ones(size(b1));

for j = 1:ncol
   ddf1(:,j) = df1;
   ddf2(:,j) = df2;
end

df1 = ddf1;
df2 = ddf2;

% since we already have standard error (se) values ...
sez = sqrt(se1.^2 + se2.^2);


dfe = (se1.^2 + se2.^2).^2 ./ (se1.^4 ./ df1 + se2.^4 ./ df2); % confirmed correct

dfer = floor(dfe); % just to be slightly more conservative

wt = (b1 - b2) ./ sez; % confirmed correct

% give the two-tailed p-value by default, using the corrected d.f.

p = 2*tcdf(-abs(wt),dfer);