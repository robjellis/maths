function out = ci_sim(P,n,width,iter)

%
% out = ci_sim(P,n,width,iter)
% 
% note: only for SINGLE PROPORTION results (e.g., proportion of successes
% ranging from 0 to 1)
%
% width = width of the CI (e.g., 80, 90, 95, 99)
%
% all methods tested have a manual correction for x == 0 or x == n
%

% RJ Ellis | 2017.11.18

% safety check
if P > 1
    error('P must be between 0 and 1, inclusive')
end

if nargin < 3
    width = 95;
end

if nargin < 4
    iter = 10000;
end

% declare variables

% does the CI *include* the true proportion? 1 = yes, 0 = no
inc_wald_nc     = zeros(iter,1);  
inc_wald_cc     = zeros(iter,1);
inc_wilson_nc   = zeros(iter,1);
inc_wilson_cc   = zeros(iter,1);

% is the *lower bound* of the CI <= the true proportion? 1 = yes, 0 = no
% This is a slightly different criteria, but relevant for the "star
% ratings" project; since we want to know if the CI is conservative, we are
% OK with lower bounds that fall below the true proportion

lb_wald_nc     = zeros(iter,1);  
lb_wald_cc     = zeros(iter,1);
lb_wilson_nc   = zeros(iter,1);
lb_wilson_cc   = zeros(iter,1);
    
% *******************************
% THE LOOP
for i = 1:iter
    
    % get the sample
    samp = rand(n,1) < P; % very simple inverse transform; 
                           % values are from open interval (0,1)
                           
    x = sum(samp); % Vollset 1993 convention

    B = bino_ci_calc(x,n,width);
    
    % store the results
    
    % need to use inequality because of potential for manual corrections if x == 0 or x == n
    
    % does CI include the true proportion?
    if B.wald_nc(1) <= P && P <= B.wald_nc(2) 
        inc_wald_nc(i) = 1;
    end
    
    if B.wald_cc(1) <= P && P <= B.wald_cc(2)
        inc_wald_cc(i) = 1;
    end

    if B.wils_nc(1) <= P && P <= B.wils_nc(2)
        inc_wilson_nc(i) = 1;
    end   

    if  B.wils_cc(1) <= P && P <= B.wils_cc(2)
        inc_wilson_cc(i) = 1;
    end
    
    % is the lower bound <= the true proportion?
    if B.wald_nc(1) <= P  
        lb_wald_nc(i) = 1;
    end
    
    if B.wald_cc(1) <= P 
        lb_wald_cc(i) = 1;
    end

    if B.wils_nc(1) <= P 
        lb_wilson_nc(i) = 1;
    end   

    if  B.wils_cc(1) <= P 
        lb_wilson_cc(i) = 1;
    end    
    
    
end

if iter == 1
   % just to report out for dummy checking
   out.walc_nc = wald_nc;
   out.wald_cc = wald_cc;
   out.wils_nc = wilson_nc;
   out.wils_cc = wilson_cc;   

else
    % as a percentage
    out.inc_wald_nc = 100 * sum(inc_wald_nc) / iter;
    out.inc_wald_cc = 100 * sum(inc_wald_cc) / iter;
    out.inc_wils_nc = 100 * sum(inc_wilson_nc) / iter;
    out.inc_wils_cc = 100 * sum(inc_wilson_cc) / iter;
    
    out.lb_wald_nc = 100 * sum(lb_wald_nc) / iter;
    out.lb_wald_cc = 100 * sum(lb_wald_cc) / iter;
    out.lb_wils_nc = 100 * sum(lb_wilson_nc) / iter;
    out.lb_wils_cc = 100 * sum(lb_wilson_cc) / iter;    
end

