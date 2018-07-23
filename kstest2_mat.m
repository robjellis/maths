function output = kstest2_mat(X, Y, mnX, mnY, nX, nY)

% special code to do many 2-sample KS tests at once
% X is 1D, N rows x 1 column
% Y is 2D, with N rows x M columns

% mnX is a scalar; the mean of the original data values
% mnY is a 1 x M vector; the means of each column of Y

% nX is a scalar; the number of items originally in X
% nY is a 1 x M vector; the number of items originally in each column of Y

% check to make sure same number of rows

if size(X,2) > 1
    % bad
    error(' X can only be one-dimensional ')
end

if size(X,1) == size(Y,1)
    % OK
else
    error(' Data structure is not correct ')
end

nrows = size(X,1);
ncols = size(Y,2);

% get X to be the same size as Y
X = repmat(X,[1 ncols]);

% same for nX
nX = repmat(nX,[1 ncols]);

% same for mnX
mnX = repmat(mnX,[1 ncols]);

% subtract Y - X; item-wise subtraction
diff = abs(Y - X);

% same for means
mn_diff = mnY - mnX;

% we only care about the direction
mn_sig = sign(mn_diff) / 2; % so values are either +.5 or -.5

% get the KS statistic for each column
D = max(diff,[],1);


% via Matlab kstest2: Compute the asymptotic P-value approximation and accept or
% reject the null hypothesis on the basis of the P-value.

n      =  nX .* nY /(nX + nY);
lambda =  max((sqrt(n) + 0.12 + 0.11./sqrt(n)) .* D, 0);

% via Matlab kstest2: one-tailed P-value; one value for each column
p_one  =  exp(-2 * lambda .* lambda);

% RJE: now we just double this and do a binary significance test
p_two = p_one * 2;

% binary significance
p_sig = p_two < .05;

% now bring in the sign of the mean
sig_adj = p_sig .* mn_sig; % will stay at 0 unless the test is significant!

% now we get the final result

sig_final = 50 + 100 * sum(sig_adj)/ncols;

output.sig_final = sig_final;





