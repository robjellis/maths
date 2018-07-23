function out = beta_gof_old(X,do_plot)

% goodness-of-fit test - THIS METHOD IS NON-OPTIMAL
% it never returns an RMSE = 0, because there is always subtle error
% between ECDF values and CDF values. Instead, better to just do a brute
% force simulation and compare Pr(Model > Data)
%
% INSTEAD USE beta_gof_sim.m
%
% RJE | 9 Feb 2018

if nargin < 2
	do_plot = 1;
end

% make sure no 0s or 1s
if sum(X==0) > 0 || sum(X==1) > 0
	fprintf(' (Removing 0s and 1s from data.)')
	X = X(and(X>0,X<1));
end

% --------------------------
% get ECDF of data
[f, X_u] = ecdf(X);

% fit a beta distribution to this data
bparam = betafit(X);
a0 = bparam(1);
b0 = bparam(2);

% get CDF values for modeled beta at original X locations
beta_cdf = betacdf(X_u,a0,b0);

% get the difference
D = beta_cdf - f;

% median abs difference
mad = median(abs(D));

% find residuals of OLS REGRESSION
% strongest test: there should be a linear fit between the CDF values
XX = [ones(size(X_u)) X_u]; % must include a dummy column of ones

[~, ~, r, ~, stats] = regress(beta_cdf,XX);
ols_R2 = stats(1);
ols_rmse = std(r,0); % flag = 0 for correction by N - 1

pears = corr(X_u,beta_cdf);
rho = corr(X_u,beta_cdf,'type','Spearman');

%% plots
if do_plot == 1
	
	figure(20)
	plot(X_u, f,'r','LineWidth',1.2)
	hold on
	plot(X_u,beta_cdf,'g')
	hold off
	xlim([0 1])
	title('ECDF of data vs CDF of Beta')
	xlabel('Value')
	ylabel('Cumulative prob.')

	figure(21)
	plot(f,beta_cdf,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1],...
		'MarkerSize',10,...
		'Marker','.',...
		'Color',[1 0 0])
	xlabel('Data EDF values')
	ylabel('Beta CDF values')
	hold on
	plot([0 1],[0 1],'g')
	hold off
	
	figure(22)
	ecdf(D)
	xlabel('Beta CDF - Data ECDF values')
end

%% outputs
out.num_unique = numel(unique(X));
out.a0 = a0;
out.b0 = b0;
out.diff_mean = mean(D);
out.diff_std = std(D);
out.diff_mad = mad;
out.ols_rmse = ols_rmse;
out.ols_R2	 = ols_R2;
out.pearson  = pears;
out.spearman = rho; % will ALWAYS be 1.0 since ECDF and CDF are monotonically related 



