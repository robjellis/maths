function [output] = corr_diff(N,iter)

% compares Z values from the Rosenthal test for correlations, versus a
% simple subtraction measure

% variables
r_A_B1 = zeros(iter,1);
r_A_B2 = zeros(iter,1);
r_B1_B2 = zeros(iter,1); 
r_A_Bd = zeros(iter,1);
diff_z = zeros(iter,1);
rosen_z = zeros(iter,1);

progressbar

for i = 1:iter
    
    progressbar(i/iter) % Update progress bar
    % random vectors
    A = randn(N,1);
    B1 = randn(N,1);
    B2 = randn(N,1);
    Bd = B1 - B2; % will also be normally distributed
    
    % actual correlations
    r_A_B1(i) = corr(A,B1);
    r_A_B2(i) = corr(A,B2);
    r_B1_B2(i) = corr(B1,B2);
    r_A_Bd(i) = corr(A,Bd);
    
    % Fisher z
    %diff_z(i) = r2z(r_A_Bd(i));
    
    % actual Z value (i.e., convert r to t, t to p, and p to z)
    [p t z] = r2p(r_A_Bd(i),'b',N,1); % note: tails doesn't matter for t-value
    diff_z(i) = z;
    
    % do the Rosenthal calculation
    [Z P] = rosenthal_Z(r_A_B1(i),r_A_B2(i),r_B1_B2(i),N);
    rosen_z(i) = Z;
end

B1_min_B2 = r_A_B1 - r_A_B2;

% plots

figure(11)
plot(B1_min_B2/2,r_A_Bd,'.');  % ** better match when we divide the difference by 2

figure(12)
hist(r_A_Bd,30);

figure(13)
hist(rosen_z,100)
xlabel('Rosenthal Z-value')

figure(14)
plot(diff_z,rosen_z,'.')
xlabel('Difference z');
ylabel('Rosenthal Z');


% difference of p-values
z_diff = rosen_z - diff_z;

% k-s tests
[rh rp] = kstest(rosen_z);
[dh dp] = kstest(diff_z);

figure(15)
hist(z_diff,30)
xlabel('Rozen Z minus diff z')

% 
% output

zcrit = norminv(.95); % one-tailed, p = .05

output.rosenthal_ks = rp;
output.rosenthal_alpha = sum(rosen_z > zcrit)/iter;
output.difference_ks = dp;
output.difference_alpha = sum(diff_z > zcrit)/iter;





    
    