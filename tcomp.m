function [stats] = tcomp(N,alpha,tail,iter)

% t-test comparison: 
% 1. paired
% 2. 2-sample, equal variance
% 3. 2-sample, unequal variance

% redefine tail as one-tailed (X - Y) or 2-tailed
if tail == 1
   tail = 'right';
elseif tail == 2
   tail = 'both';
end

% set up variables


for i = 1:iter  
    x = randn(N,1);
    y = randn(N,1);
    
    % paired t-test
    [a b c d] = ttest(x,y,alpha,tail);
    tP_sig1(i) = a;
    tP_stat(i) = d.tstat;
    
    % 2-sample, equal variance
    [a b c d] = ttest2(x,y,alpha,tail,'equal');
    tE_sig1(i) = a;
    tE_stat(i) = d.tstat;
    
    % 2-sample, unequal variance
    [a b c d] = ttest2(x,y,alpha,tail,'unequal');
    tU_sig1(i) = a; % 1 if result is significant
    tU_stat(i) = d.tstat; % the actual T-statistic
    
end

% now, we want to see, given the *same actual T-crit*, which survives

tcrit = tinv(1-alpha,N-1);

tP_sig2 = sum(tP_stat >= tcrit);
tE_sig2 = sum(tP_stat >= tcrit);
tU_sig2 = sum(tP_stat >= tcrit);
