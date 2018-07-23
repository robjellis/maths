function [output] = sl_corr(N,iter)

% signed log effect on correlation

% variables
r_A_B = zeros(iter,1);
r_slA_slB = zeros(iter,1);

p_A_B = zeros(iter,1);
p_slA_slB = zeros(iter,1);

for i = 1:iter
    
    % get the vectors
    %A = randn(N,1);
    %B = randn(N,1);

    A = exprnd(1,N,1) - 1; % mean will be 0
    B = exprnd(1,N,1) - 1; % mean will be 0
    
    % transform them
    slA = sign(A) .* log(abs(A)+1);
    slB = sign(B) .* log(abs(B)+1);

    % correlations
    [r p] = corr(A,B);
    r_A_B(i) = r;
    p_A_B(i) = p;

    [r p] = corr(slA,slB);
    r_slA_slB(i) = r;
    p_slA_slB(i) = p;

end

output.perc_p_A_B = sum(p_A_B < .05)/iter*100;
output.perc_p_slA_slB = sum(p_slA_slB < .05)/iter*100;
