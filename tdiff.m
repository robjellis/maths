function [res] = tdiff(N1,N2,iter,pcrit)

% compare "tdiff" values (independent samples) versus the result of a
% traditional two-sample t-test

% create the theoretical diff distribution
x1 = trnd(N1-1,10000000,1);
x2 = trnd(N2-1,10000000,1);
zdiff = t2z(x1,N1-1) - t2z(x2,N2-1);

%x1 = randn(1000000,1);
%x2 = randn(1000000,1);
%zdiff = x1 - x2;

diff_crit = prctile(zdiff,100*(1-pcrit));

% and the critical two-sample value

two_crit = t2z(tinv(1 - pcrit, N1+N2-2),(N1+N2-2));
    
    r1 = randn(N1,iter);
    r2 = randn(N2,iter);
    
    % two one-sample tests
    [H,P,CI,stats] = ttest(r1);
    t1 = stats.tstat;
    
    [H,P,CI,stats] = ttest(r2);
    t2 = stats.tstat;
    
    % convert t to z so that things aren't as extreme!
    z1 = t2z(t1,N1-1);
    z2 = t2z(t2,N2-1);
    
    z_diff = z1 - z2;
    
    % two-sample test
    
    [H,P,CI,stats] = ttest(r1,r2);
    t_two = stats.tstat;
    
    z_two = t2z(t_two,(N1+N2-2));
    
% significant? (1 = yes)
z_two_sig = z_two >= two_crit;
z_diff_sig = z_diff >= diff_crit;

% what we really want is a contingency table

cor_reg  = (z_two_sig == 0 & z_diff_sig == 0);
hit      = (z_two_sig == 1 & z_diff_sig == 1);
miss     = (z_two_sig == 1 & z_diff_sig == 0);
false_al = (z_two_sig == 0 & z_diff_sig == 1);

% now we divide up the voxels by color

figure(50)
plot(z_two(cor_reg == 1),z_diff(cor_reg == 1),'g+')
hold on
plot(z_two(hit == 1),z_diff(hit == 1),'b+')
plot(z_two(miss == 1),z_diff(miss == 1),'m+')
plot(z_two(false_al == 1),z_diff(false_al == 1),'r+')
hold off


res.cor_reg  = sum(cor_reg)/iter;
res.hit      = sum(hit)/iter;
res.miss     = sum(miss)/iter;
res.false_al = sum(false_al)/iter;
res.HplusFA  = res.hit+res.false_al;

% sum(t_diff_sig)/iter
% dprime = norminv(hit) - norminv(false_al)


    
    