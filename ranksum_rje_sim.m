function out = ranksum_rje_sim(n,iter)

% how close to .05 do we get?

H_rje = zeros(iter,1);
H_mat = zeros(iter,1);
H_ttest = zeros(iter,1);

for i = 1:iter
    x = randn(n,1);
    y = randn(n,1); % independent
    
    do_exact = 0;
    rje = ranksum_rje(x,y,do_exact);
    
    H_rje(i) = rje.H;
    
    [p h] = ranksum(x,y);
    H_mat(i) = h;
    
    H_ttest(i) = ttest2(x,y);
    
end

out.FPR_rje = sum(H_rje)/iter;
out.FPR_mat = sum(H_mat)/iter;
out.FPR_ttest = sum(H_ttest)/iter;



