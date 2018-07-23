function [supra tA tB tminAB tN] = null_ttests(nsub,nvox)

%
% function [tA tB tminAB tN] = null_ttests(nsub,nvox)
%
% pvox = voxel-level p-value (e.g., p = .005)
% default: ncond = 2 for now

% 1. get the arrays for condA and condB

condA = randn(nsub,nvox);
condB = randn(nsub,nvox);

meanAB = (condA + condB) / 2;

condN = nullit(condA,condB);

%figure
%hist(condN(:),50)

% 2. actual t-ttests

[w x y z] = ttest(condA);
tA = z.tstat;

    minA = min(tA);
    maxA = max(tA);

[w x y z] = ttest(condB);
tB = z.tstat;

    minB = min(tB);
    maxB = max(tB);

[w x y z] = ttest(meanAB);
tmeanAB = z.tstat;

[w x y z] = ttest(condN);
tN = z.tstat;

% 3a. min statistic at 2nd-level

tminAB = min(tA,tB);

%figure
%plot(tminAB(:),tN(:),'.')

% how many supra-threshold voxels are there? at pvox

pvox = [.5 .1 .05 .01 .005 .001 .0005];

tcrit = abs(tinv(pvox,11));

for t = 1:numel(tcrit)
   
    sig_tA(t) = sum(tA>tcrit(t));
    sig_tB(t) = sum(tB>tcrit(t));
    sig_tminAB(t) = sum(tminAB>tcrit(t));
    sig_tnullAB(t) = sum(tN>tcrit(t));
end

supra.pvals = pvox;
supra.sig_tA = sig_tA;
supra.sig_tB = sig_tB;
supra.sig_tminAB = sig_tminAB;
supra.sig_tnullAB = sig_tnullAB;

suprathreshold_voxels = supra


% 4. histograms

%minplot = min(minA,minB);
%maxplot = max(maxA,maxB);

minx = -6;
maxx = 6;

xvals = linspace(minx,maxx,50);

figure(51)
hist(tA,xvals)
title('cond_A, t-tests')
axis([minx maxx 0 20000])

figure(52)
hist(tB,xvals)
title('cond_B, t-tests')
axis([minx maxx 0 20000])


%figure(53)
%hist(tmeanAB,xvals)
%title('avgerage of cond_A and cond_B, t-tests')
%axis([minx maxx 0 20000])

figure(53)
hist(tminAB,xvals)
title('second-level minimum statitic, t-values')
axis([minx maxx 0 20000])

figure(54)
hist(tN,xvals)
title('first-level null statistic, t-tests')
axis([minx maxx 0 20000])
