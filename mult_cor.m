function [stats] = mult_cor(alpha,svals,pvals)

% function [stats] = mult_cor(alpha,svals,pvals)
%
% various multiple correction methods
% alpha = the desired FWE p or FDR q value
% svals = the actual statistic values (t-values, z-values, binomial); can
%         be a 3D volume
% pvals = the associated p-values; can be a 3D volume
%
% The program will *adjust* the actual P-values to their "corrected" level,
% per the different formulas.  Any P_adj that are <= alpha are then deemed
% significant.
%
% for a full description of the FWE tests and their calculations, see
% http://support.sas.com/documentation/cdl/en/statug/63347/HTML/default/viewer.htm#statug_multtest_sect014.htm
%
% all output has been validated against the "R" stats package: 
% http://stat.ethz.ch/R-manual/R-patched/library/stats/html/p.adjust.html
%
% last update = 2012.01.13
%
% Copyright Robert J Ellis

svals2 = svals(:);
pvals2 = pvals(:);

svals = svals2(isfinite(pvals));  % Toss NaN's (e.g., non-brain voxels)
pvals = pvals2(isfinite(pvals));

nums = numel(svals);
nump = numel(pvals);

% of course, we need nums and nump to be equal
if (nums - nump) ~= 0
   fprintf('\n Warning: number of statistic values and p-values is not equal. Terminating.\n\n');
   return
end

% first, sort the s-vals and p-vals together

sp = [svals pvals];
sp_asc = sortrows(sp, 2);   % sort by p-vals, ascending

s = sp_asc(:,1);            % pull out the sorted s and p as vectors        
p = sp_asc(:,2);

% ------------
% 1. traditional Bonferroni

corp = p * nump;  % inflate all p-vals
corp(corp>1) = 1;

Bon_p = p(max(find(corp <= alpha)));             % the smallest observed sig. P-val (can be empty)


if isempty(Bon_p)
   Bon_p = NaN;
   Bon_s = Inf;                       % no valid threshold!
else
   Bon_s = min(s(p <= Bon_p)); 
end

% ------------
% 2. Bonferroni-Holm: "step-down" procedure

mult = nump:-1:1;                    % to multiply p-observed values
mult = mult(:);
comp = mult .* p; 

corp = zeros(nump,1);
for i = 1:nump
    if i == 1
       corp(i) = comp(i);
    else
       corp(i) = max(corp(i-1),comp(i)); 
    end
    
end

corp(corp>1) = 1;            % any adjusted P larger than 1 is set to 1

BonHolm_p  = p(max(find(corp <= alpha))); % the observed P-value cutoff

BonHolmP = corp;  % all corrected Pvals

if isempty(BonHolm_p)                 % there is no valid correction
    BonHolm_p = NaN; 
    BonHolm_s = Inf;
else                                  % we are OK   
    BonHolm_s = min(s(p <= BonHolm_p));          % the smallest significant s-value
end

% ------------
% 3. Bonferroni-Hochberg: "step-up" procedure

mult = 1:nump;                       % multiplier, ascending
mult = mult(:);
corp = mult .* p;                       

for i = nump:-1:1
    if i == nump
       corp(i) = comp(i);
    else
       corp(i) = min(corp(i+1),comp(i)); 
    end
end

corp(corp>1) = 1;

BonHoch_p  = p(max(find(corp <= alpha))); % the observed P-value cutoff

BonHochP = corp;   % all corrected Pvals
  
if isempty(BonHoch_p)                 % there is no valid correction
    BonHoch_p = NaN; 
    BonHoch_s = Inf;
else                                             % we are OK   
    BonHoch_s = min(s(p <= BonHoch_p));          % the smallest significant s-value
end

% comparison of Bon-Holm and Bon-Hoch procedures

%diff_meth = BonHolmP - BonHochP;
%figure(101)
%plot(p,diff_meth,'.');

%diff_meth_mode = mode(diff_meth);
%diff_meth_max = max(diff_meth);

diff_meth_corr = corr(BonHolmP,BonHochP);


% ------------endsB_vol_Bon > 0
% 4. Benjamini-Hochberg FDR

q = alpha;

% the following is by Tom Nichols
% http://www2.warwick.ac.uk/fac/sci/statistics/staff/academic-research/nichols/software/fdr

    %{

    This function takethr_FDR_all = zeros(numf,1);
    thr_Prop_all = zeros(numf,1);s a vector of p-values and a FDR rate. 
    It returns two p-value thresholds, one based on an assumption of independence or positive dependence, 
    and one that makes no assumptions about how the tests are correlated. For imaging data, an assumption of 
    positive dependence is reasonable, so it should be OK to use the first (more sensitive) threshold.

    %}

    % function [pID,pN] = FDR(p,q)
    % FORMAT [pID,pN] = FDR(p,q)
    % 
    % p   - vector of p-values
    % q   - False Discovery Rate level
    %
    % pID - p-value threshold based on independence or positive dependence
    % pN  - Nonparametric p-value threshold

    % p = p(fininte(p));  % Toss NaN's   
    % p = sort(p);        % p is already sorted above
    
    V = length(p);
    I = (1:V)';

    cVID = 1;
    cVN = sum(1./(1:V));
    
    corp = I/V*q/cVID;                % the FDR critical value at each test point
    
    FDR_p = p(max(find(p <= corp)));  % calculation for independent or pos. correlations

    
    
if isempty(FDR_p)
   FDR_p = NaN;                           % no observed values reach FDR significance!
   FDR_s = Inf;
else
    % we are OK
    FDR_s = min(s(p <= FDR_p));
end


% ------------
% get all the results

stats.Bon_p      = Bon_p;
stats.Bon_s      = Bon_s;
stats.BonHolm_p  = BonHolm_p;
stats.BonHolm_s  = BonHolm_s;
stats.BonHoch_p  = BonHoch_p;
stats.BonHoch_s  = BonHoch_s;
stats.diff_meth_corr = diff_meth_corr;
stats.FDR_p      = FDR_p;
stats.FDR_s      = FDR_s;

