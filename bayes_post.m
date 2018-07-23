function post = bayes_post(A_suc,A_tot,B_suc,B_tot)

% THIS IS A FLAWED METHOD; WAS JUST CURIOUS TO SEE WHAT IT WOULD DO!!

% calculate posterior probability p(B|suc) given:
%   A_suc = number of successes in A
%   A_tot = total number of A events
%   B_suc = number of successes in B
%   B_tot = total number of B events
%
% RJE | April 2018

prob_B           = B_tot / (A_tot + B_tot);
prob_suc_given_B = B_suc / B_tot;
prob_suc         = (A_suc + B_suc) / (A_tot + B_tot);

post = (prob_B * prob_suc_given_B) / prob_suc;