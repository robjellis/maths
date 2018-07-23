efunction out = ranksum_1samp_sim(N,add,iter,meth)

% store
sig_mw  = zeros(iter,1);
sig_t1  = zeros(iter,1);
sig_t2  = zeros(iter,1);

progressbar(0)

% determine bounds for simulated
lv = (1 / N) / 2;
uv = 1 - lv;

if strcmp(meth,'quant')
    % we have a single sample 2
    samp2 = norminv(linspace(lv,uv,N)); % can't use 0 and 1 as bounds since those yield -Inf and +Inf, respectively
end

for i = 1:iter
   
   samp1 = randn(N,1) + add;
   
   if strcmp(meth,'randn')
       % new vector each time
       samp2 = randn(N,1);
   end
   
   % do the test
   [p, h] = ranksum(samp1,samp2);
   
   sig_mw(i) = h; % 0 or 1
   
   % one-sample t
   sig_t1(i) = ttest(samp1);
   
   % two-sample t
   sig_t2(i) = ttest2(samp1,samp2);
   
   progressbar(i/iter)
end

progressbar(1)

% output
out.sig_mw_prc  = 100 * sum(sig_mw) / iter;
out.sig_t1_prc  = 100 * sum(sig_t1) / iter;
out.sig_t2_prc  = 100 * sum(sig_t2) / iter;
if strcmp(meth,'quant')
    out.samp2       = samp2;
end