function out = corr_subsample(X,Y,N,iter)

% assume X and Y are a "population" of values that we take random subsamples "x"
% and "y" from, where numel(x) = numel(y) = N
%
% Thee purpose: to see how the Spearman rho values of the subsamples
% compare with the population rho 

format short
% population rho
pop_rho = corr(X,Y,'type','Spearman');

nX = numel(X);

rho = nan(nX,1);
mad = nan(nX,1);

progressbar(0)
for i = 1:iter
	% random permutation
	ind = randperm(nX);
	
	% take first
	ind = ind(1:N);
	
	% get the data
	x = X(ind);
	y = Y(ind);
	
	% find their ranks
	rX = tiedrank(x);
	rY = tiedrank(y);
	
	% take the abs diff
	d = abs(rY - rX);
	
	% what is the median absolute change in rank?
	mad(i) = median(d);
	
	% correlation
	rho(i) = corr(x,y,'type','Spearman');
	
	if rem(i,10) == 0
		progressbar(i/iter)
	end
	
end

progressbar(1)

figure(20)
ecdf(rho)
xlabel('Spearman rho')

figure(21)
ecdf(mad)
xlabel('Median absolute change in rank')

out.pop_rho = pop_rho;