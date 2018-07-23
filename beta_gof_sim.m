function out = beta_gof_sim(X,a0,b0,simN,do_plot)

% goodness-of-fit test 
% brute force simulation to evaluate compare Pr(Model > Data)
%
% saves time by simply computing the KS 2-sample and Mann-Whitney PS statistics on simulated
% values versus real values in X, thus not requiring iterations. For
% iteration-based version, see beta_gof_sim_old.m
%
% RJE | 9 Feb 2018

tic;

X = X(:);

% min and max values of data
X_min = min(X);
X_med = median(X);
X_max = max(X);

if nargin < 3
	calc_beta = 1;
else
	calc_beta = 0;
end

if nargin < 4
	simN = 1000000;
end

if nargin < 5
	do_plot = 1;
end

% make sure no 0s or 1s
if sum(X==0) > 0 || sum(X==1) > 0
	fprintf(' (Removing 0s and 1s from data before modeling.)')
	X = X(and(X>0,X<1));
end

% fit a beta distribution to this data
if calc_beta == 1
	[bparam, pci] = betafit(X);
	a0 = bparam(1);
	b0 = bparam(2);
else
	% we input a0 and b0
end

% simulate values from new distribution
bvals = betarnd(a0,b0,simN,1);
b_min = min(bvals);
b_med = median(bvals);
b_max = max(bvals);

% different analysis
[ff, xx]  = ecdf(X); % get ECDF of data	
beta_cdf2 = betacdf(xx,a0,b0); % CDF values matched to these data values

% diffference
D = beta_cdf2 - ff;

% med abs diff
med_AD = median(abs(D));
max_AD = max(abs(D));

% Mann Whitney PS value
mw = ranksum_rje(X,bvals);

% KS test 2-sample
[~, KS_pval, KS_stat] = kstest2(X,bvals);

% AD test 2-sample
% http://www.jaqm.ro/issues/volume-6,issue-3/pdfs/1_engmann_cousineau.pdf
% Github: https://github.com/vncntprvst/tools/blob/master/adtest2.m
% Problem: WAY TOO SLOW for large samples like we have here

% note: code is written to accept row vectors
%[H, adstat, critvalue] = adtest2(X',bvals')

tocc = toc;

%% plot

if do_plot == 1
	% full range of x-values
	xfull = [0:.01:1]'; % column
	

	
	% ECDF of simulated beta
	[fb, xb] = ecdf(bvals);
	
	% we actually want to show the ECDF for the full range of 0 to 1
	xx_lo = xfull(xfull < X_min);
	xx_hi = xfull(xfull > X_max);
	
	ff_lo = zeros(size(xx_lo));
	ff_hi = ones(size(xx_hi));
	
	xx_mod = [xx_lo; xx; xx_hi]; 
	ff_mod = [ff_lo; ff; ff_hi];
	
	% get PDF and CDF values for modeled beta
	beta_pdf = betapdf(xfull,a0,b0);
	beta_cdf = betacdf(xfull,a0,b0);
	
	thr = .0001; % to make the plots nicer

	figure(30)
	plot(xfull(beta_pdf>thr),beta_cdf(beta_pdf>thr),'g')
	hold on
	plot(xx_mod, ff_mod,'r','LineWidth',1.2)
	plot(xb,fb,'b')
	
	hold off
	xlim([0 1])
	title('ECDF of data vs CDF of Beta')
	xlabel('Value')
	ylabel('Cumulative prob.')
	minx = min(xfull(beta_pdf>thr));
	maxx = max(xfull(beta_pdf>thr));
	xlim([minx maxx])
	

	figure(31)
	plot(ff,beta_cdf2,'.')

	corr(ff,beta_cdf2)

end

%% outputs
format short
out.a0			= a0;
out.b0			= b0;

out.www			= '[min, med, max]';
out.data		= [X_min X_med X_max];
out.sim_beta	= [b_min b_med b_max];

out.xxx			= 'KS 2-sample test';
out.med_AD      = med_AD;
out.max_AD      = max_AD;
out.KS_stat		= KS_stat;
out.KS_pval		= KS_pval;


out.yyy			= 'Pr(Model > Data)';
out.MW_PS		= mw.PS;

out.duration_sec = tocc;



