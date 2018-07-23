function out = quant_beta(X,P,Qobs)

% RJE method of quantile matching based search for optimal a0 and b0
% parameters from a Beta distribution matched to data in X, or directly to
% observed quantiles in Qobs

tic;

if nargin < 2
	P = [.01 .50 .99];
	%P = [.01 .10 .25 .50 .75 .90 .99]; % the target CDF values to match on
end

if nargin < 3
	X = X(:);
	Qobs = quantile(X,P);
end

% determine overall shape
minQ = min(Qobs);
medQ = Qobs(P == 0.5);
maxQ = max(Qobs);

dH = maxQ - medQ;
dL = medQ - minQ;

% log ratio
lr = log2(dH/dL);

if lr > 0.5
	skew = 1; % right skew (positive); only test b0 values that are >= a0
elseif lr < -0.5
	skew = -1; % left skew (negative); only test b0 values that are <= a0
else
	skew = 0; % minimal skew
end

skew

% default max likelihood
param = betafit(X);

% 0. set up

numP = numel(P);

% candidate a0 and b0 values
minV = .5;
maxV = 100;
step = .5;

% define a0 up front
a0 = minV:step:maxV;

tot = numel(a0)^2; % total possible

% 1. get observed quantiles from data in X



% 2. large table
a0_test = nan(tot,1);
b0_test = nan(tot,1);

maxAE = nan(tot,1);

ctr = 1;

for i = 1:numel(a0)
	
	% get b0 values to test
	if skew == -1
		b0 = minV:step:a0(i);
	elseif skew == 1
		b0 = a0(i):step:maxV;
	elseif skew == 0
		b0 = minV:step:maxV;
	end
	
	for j = 1:numel(b0)
		a0_test(ctr) = a0(i);
		b0_test(ctr) = b0(j);
		
		Qtheo = betainv(P,a0(i),b0(j));
		
		maxAE(ctr) = max(abs(Qobs - Qtheo));
		
% 		if MAE(ctr) < .02
% 			% only save results that improve the fit
% 			ctr = ctr + 1;
% 		else
% 			% don't save
% 		end
		
		ctr = ctr + 1;
	end
end

% remove empty


T = table(a0_test,b0_test,maxAE);
T = T(1:ctr-1,:);

ind = maxAE == min(maxAE);

a0_best = T.a0_test(ind);
b0_best = T.b0_test(ind);

tocc = toc;

%% figures

if numel(X) > 0
	figure(400)
	ecdf(X)
	xlim([0 1])
	xlabel('Data values')
else
	figure(400)
	close
end

% figure(410)
% plot(maxAE)
% ylabel('Mean Abs Diff')

figure(420)
pointsize = 50;
scatter(T.a0_test,T.b0_test,pointsize,-log10(T.maxAE),'.')


%% outputs
out.P_eval	= P;
out.Qobs	= Qobs;
out.T		= T;
out.quant_a0 = a0_best;
out.quant_b0 = b0_best;
out.quant_maxAE = min(maxAE);
out.ml_a0	= param(1);
out.ml_bo	= param(2);
out.duration_sec = tocc;
