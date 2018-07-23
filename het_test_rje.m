function het_test_rje(x,y)

% testing
% note: X is assumed to be a single vector for simplicity

% HetTest.m won't work for very large datasets!

if numel(x) > 15000
	warning(' Underlying code (HetTest.m) struggles with large vectors of data.\n Suggest taking (multiple) subsamples instead.')
end

% 1. Do the regression
xx = [ones(size(x)) x];
[b, bint, r, rint, stats] = regress(y,xx);

% 2. plots
figure(50)
subplot(1,3,1)
plot(x,y,'.')
xlabel('X')
ylabel('Y')
title(['X vs Y: r = ' num2str(corr(x,y))])

subplot(1,3,2)
plot(x,r,'.')
xlabel('X')
ylabel('Resid(Y)')
title(['X vs Resid(Y): r = ' num2str(corr(x,r))])

subplot(1,3,3)
plot(y,r,'.')
xlabel('Y')
ylabel('Resid(Y)')
title(['Y vs Resid(Y): r = ' num2str(corr(y,r))])

% 3. test for hetero
% PVAL = TestHet(RES, X, WHICHTEST, YHAT)
het_pval = TestHet(r,x,'-BPK')
