function out = beta_gof_sim_old(X,iter,do_plot)

% goodness-of-fit test 
% brute force simulation to evaluate compare Pr(Model > Data)
%
% RJE | 9 Feb 2018

beep off

tic;

if nargin < 2
	iter = 5000;
end

if nargin < 3
	do_plot = 1;
end

% make sure no 0s or 1s
if sum(X==0) > 0 || sum(X==1) > 0
	fprintf(' (Removing 0s and 1s from data before modeling.)')
	X = X(and(X>0,X<1));
end

nX = numel(X);

X_u = unique(X);

Pr = nan(iter,1);

% fit a beta distribution to this data
bparam = betafit(X);
a0 = bparam(1);
b0 = bparam(2);

progressbar(0)

for i = 1:iter

	% simulate values from new distribution
	bvals = betarnd(a0,b0,nX,1);

	% compute rate of Beta > Data
	d = bvals - X;

	Pr(i) = mean(d>0); % taking average is same as taking sum and dividing by iter
	
	
	% alternative - don't use this	
% 	d2 = sort(bvals) - sort(X);
% 	
% 	std_d(i) = std(d2);
% 	
% 	figure(14)
% 	plot(sort(X),sort(bvals),'.')
% 	
% 	figure(15)
% 	ecdf(d2)

	if rem(i,100) == 0
		progressbar(i/iter)
	end
end

progressbar(1)

% get the interval
Pr_lo = prctile(Pr,2.5);
Pr_med = prctile(Pr, 50);
Pr_hi = prctile(Pr,97.5);
tocc = toc;

%% plot

if do_plot
	% get ECDF of data
	[ff, xx] = ecdf(X);

	% get CDF values for modeled beta
	beta_cdf = betacdf(X_u,a0,b0);
	
	thr = .0001;

	figure(30)
	plot(xx(ff>thr), ff(ff>thr),'r','LineWidth',1.2)
	hold on
	plot(X_u(beta_cdf>thr),beta_cdf(beta_cdf>thr),'g')
	hold off
	xlim([0 1])
	title('ECDF of data vs CDF of Beta')
	xlabel('Data')
	ylabel('Beta')
	minx = min([min(xx(ff>thr)) min(X_u(beta_cdf>thr))]);
	maxx = max([max(xx(ff>thr)) max(X_u(beta_cdf>thr))]);
	xlim([minx maxx])

	figure(32)
	ecdf(Pr)
	xlabel('Pr(Model > Data)')
	xlim([.4 .6])
	
end
out.num_data = nX;
out.a0 = a0;
out.b0 = b0;
out.xxx = 'Pr(Model > Data)';
out.Pr_lo = Pr_lo;
out.Pr_med = Pr_med;
out.Pr_hi = Pr_hi;
out.duration_sec = tocc;


