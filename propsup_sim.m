function out = propsup_sim(n1range,n2,dist,iter)
% out = propsup_sim([n1min n1max],n2,dist,iter)
%
% dist is either 'u' or 'n'
% let's see how distribution point estimates compare with PS values

PS1 = nan(iter,1);
Z1  = nan(iter,1);
N1  = nan(iter,1);
P   = nan(iter,1);
MD1 = nan(iter,1);
MN1 = nan(iter,1);
MN2 = nan(iter,1);

progressbar(0)
for i = 1:iter
    
    % samp1 will be the sample that varies more
    n1 = randint(1,1,[n1range(1) n1range(2)]);
    
    if strcmp(dist,'usig') % some signal is present

        c = rand(1,1); % random from uniform dist
        d = rand(1,1);

        a = min(c,d);
        b = max(c,d);
        
        samp1 = a + (b-a).*rand(n1,1);
        samp2 = rand(n2,1);
        label = 'Uniform (with signal)';
        
    elseif strcmp(dist,'unull')
        samp1 = rand(n1,1);
        samp2 = rand(n2,1);
        label = 'Uniform (null)';
        
    elseif strcmp(dist,'n')
        
        %samp1 = randn(n,1) + randn;
        samp1 = randn(n1,1);
        samp2 = randn(n2,1);
        label = 'Normal (null)';
    end
    
    N1(i) = n1;
    MN1(i) = mean(samp1);
    MD1(i) = median(samp1);

    MN2(i) = mean(samp2);
        
    res = ranksum_rje(samp1,samp2,'method','approximate');
    PS1(i) = res.PSsamp1;
    Z1(i)  = res.Zsamp1_ranksum;
    P(i)   = res.Pval_ranksum;
    
    
    % let's calculate Cohen's d
    progressbar(i/iter)
end

% let's calculate the TRUE "critical value" of PS that corresponds with a 
% Z = 1.96 (i.e., two-tailed P < .05)

zcrit1 = 1.645;  % one tailed sig; i.e., test of samp1 > samp2, which has a positive z value
                 % normcdf(1.645) = 0.950015094460879; i.e., it is *just* significant with p < .05
zcrit2 = 1.960; % this is two-tailed sig with p < .05

sigU  = sqrt(n1*n2*(n1+n2+1)/12);
mnU   = n1*n2/2;

Ucrit1  = zcrit1*sigU + mnU;
PScrit1 = Ucrit1/(n1*n2);

Ucrit2  = zcrit2*sigU + mnU;
PScrit2 = Ucrit2/(n1*n2);

progressbar(1)
% plots
opt = 5;

figure(101)
plot(MN1,PS1,'.');
xlabel(['MN1 ' label])
ylabel('PS1')

figure(102)
plot_eda(MN1,PS1,opt,102);

figure(103)
plot(MD1,PS1,'.');
xlabel(['MD1 ' label])
ylabel('PS1')

figure(104)
plot_eda(MD1,PS1,opt,104);

figure(105)
plot(0.5+MN1-MN2,PS1,'.')
xlabel('0.5 + (MN1 - MN2)')
ylabel('PS1')

figure(106)
plot(Z1,PS1,'.')
xlabel('Z1')
ylabel('PS1')
out = [];
% 
% figure(107)
% subplot(1,2,1)
% plot(N1,Z1,'.')
% xlabel('N1')
% ylabel('Z1')
% 
% subplot(1,2,2)
% plot(N1,PS1,'.')
% xlabel('N1')
% ylabel('PS1')

figure(109)
ecdf(PS1)
xlabel('PS1')

out.PS_theo_one_tail = PScrit1;
out.PS_theo_two_tail = [1 - PScrit2 PScrit2];
out.PS_obs_CI = [prctile(PS1,2.5) prctile(PS1,97.5)];

% let's also figure out the lowest and highest PS that we get *without* being significant
% out.PS_Pnon_lo = min(PS1(P >= .05));
% out.PS_Pnon_hi = max(PS1(P >= .05));
