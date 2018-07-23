function out = rje_2prop_ci(a,m,b,n,width)

%
%
% This *only* works to find the conservative "inner boound" (i.e., value closest to 0); can't use it to
% compute an upper bound as well
%
%

% where
% a = successes for group A
% m = total cases for group A
% b = successes for group B
% n = total cases for group B
% width: e.g., 95 for 95% CI

% compute the Wilson score - use continuity corrected version
resA = bino_ci_calc(a,m,width);
wA   = resA.wils_cc;

resB = bino_ci_calc(b,n,width);
wB   = resB.wils_cc;

pA = a/m;
pB = b/n;

% determine CI overlap
if wA(1) <= wB(1) && wA(2) >= wB(1)
    overlap = 1; 
elseif wA(1) <= wB(2) && wA(2) >= wB(2)
    overlap = 1;
elseif wA(1) >= wB(1) && wA(2) <= wB(2)
    overlap = 1;
elseif wB(1) >= wA(1) && wB(2) <= wA(2)
    overlap = 1;
else
    overlap = 0;
end

if overlap == 1
    % conservative approach: if any partial overlap between the CIs, then
    % it's possible there is no difference between them
    d_score = 0;
    r_score = 0;
else
    
   % **** difference in props *****
   d11 = wA(1) - wB(1);
   d12 = wA(1) - wB(2);
   d21 = wA(2) - wB(1);
   d22 = wA(2) - wB(2);
   
   d = [d11 d12 d21 d22];
   
   % **** log2 ratio of props *****
   r11 = log2(wA(1) / wB(1));
   r12 = log2(wA(1) / wB(2));
   r21 = log2(wA(2) / wB(1));
   r22 = log2(wA(2) / wB(2));
   
   r = [r11 r12 r21 r22];
   
   % can't do a z-test in this same fashion, as it directly takes the
   % counts as inputs, and this "adjusted" method screws that up
   
   % **** find values closest to 0 ****
   d_score = d(abs(d) == min(abs(d)));
   r_score = r(abs(r) == min(abs(r)));
    
end

out.pA        = pA;
out.pB        = pB;
out.wilson_A  = wA;
out.wilson_B  = wB;
out.xxx       = '------';
out.diff      = pA - pB;
out.diff_score   = d_score;
out.log2_ratio   = log2(pA / pB);
out.log2_score   = r_score;


