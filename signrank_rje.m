function output = signrank_rje(ex,co)

% Wilcoxon signed rank test for paired samples
%
% performs a paired, two-sided test of the hypothesis
% that the difference between the matched samples in the vectors X and Y
% comes from a distribution whose median is zero.
% 
% modified from Matlab's signrank.m
%
% RJE | 2014.05.27


% 0. make sure vectors are the same
if numel(co) == 1;
    co = zeros(size(ex,1),1) + co;
end

if size(ex,1) == size(co,1)
    % OK
else
    error(' co and ex are not the same size! ')
end

nc = size(ex,2);

%% loop

subs = nan(1,nc);
nonzero = nan(1,nc);
nontie = nan(1,nc);
Tn = nan(1,nc);
Tp = nan(1,nc);
Tpr = nan(1,nc); % rje formula
Pp = nan(1,nc);
Zc = nan(1,nc); % Conover Z value
Zs = nan(1,nc); % Siegel Z value
Zm = nan(1,nc); % per Matlab
sig = nan(1,nc);
Papprox = nan(1,nc);
Pexact = nan(1,nc); 


for c = 1:nc
    
    exc = ex(:,c);
    
    % 1. get the pair differences
    diffxy = exc(:) - co(:);

    % 2. cut out diff = 0 pairs
    nodiff = find(diffxy == 0);
    diffxy(nodiff) = [];
    n = length(diffxy); % this is used to get the critical value

    if n > 5

        % 3. get the ranks (on the reduced data)
        [tierank, tieadj] = tiedrank(abs(diffxy));

        % how many non-ties?
        notie = numel(unique(tierank));

        % 4. Get the rank sum
        % Matlab: "Compute signed rank statistic (most extreme version)"
        % rje: see http://en.wikipedia.org/wiki/Wilcoxon_signed-rank_test; this
        % version calculates the "T" statistic by Siegel, where only thesmaller
        % of the two sums of ranks of given sign are summed, which is easier to
        % calculate by hand

        neg = find(diffxy < 0);
        pos = find(diffxy > 0);
        Tneg = sum(tierank(neg));
        Tpos = sum(tierank(pos));
        
        % what proportion of non-zero pairs are positive?
        Ppos = numel(pos) / n;

        % Conover "T" value, using sum of ALL signed ranks (valid for ties as well)
        srank = tierank .* sign(diffxy);
        Z_con = sum(srank) / sqrt(sum(srank).^2);

        % Siegel Z value
        Z_sie = (Tpos - n*(n+1)/4) / sqrt(n*(n+1)*(2*n+1)/24);

        % RJE: the *real* value of W
        Wsign = sum(tierank.*sign(diffxy));
        Wabs = abs(Wsign);
        Wmin = min(Wabs, n*(n+1)/2-Wabs); % rje: note: this is identical to sum(1:n) - Wabs

        % do the Matlab test for clarity
        [pa h statsa] = signrank(exc,co,'method','approximate'); % gives the z value
        pe = signrank(exc,co,'method','exact'); % gives exact p value

        % 2. Matlab always returns W as the *minimum* of W and Wc, so we are already for the test

        % below are the TWO-TAILED critical values for alpha = .05, for 
        % n = 1 to n = 50 (n is the number of non-zero pair diffs)

        % sources:
        % http://web.anglia.ac.uk/numbers/biostatistics/wilcoxon/local_folder/critical_values.html
        % http://www.york.ac.uk/depts/maths/tables/wilcoxon.pdf
        % http://www.docstoc.com/docs/106461562/Critical-Values-of-the-Wilcoxon-Signed-Ranks-Test-Two-Tailed-Test

        % RJE: per Gibbons Table H, these values are signiificant for
        % two-tailed alpha (i.e., a P < .025); thus, we use <= when comparing to this!

        % in fact, these are the critical values for T- (i.e., abs value of sum of negative ranks); we just take n*(n+1)/2 minus T- to get the critical value for T+
        W_crit = [NaN NaN NaN NaN  0 0 2 3 5 8 10 13 17 21 25 29 34 40 46 52 58 65 73 81 89 98 107 116 126 137 147 159 170 182 195 208 221 235 249 264 279 294 310 327 343 361 378 396 415 434];

        % the test is significant if W is < W_crit
        sig = Wmin < W_crit(n);

        %% outputs
        subs(c) = numel(co);
        nonzero(c) = n;
        nontie(c) = notie;
        Tn(c) = Tneg;
        Tp(c) = Tpos;
        Tpr(c) = Tpos / sum(1:n); % rje formula
        Pp(c) = Ppos;
        Zc(c) = Z_con; % Conover Z value
        Zs(c) = Z_sie; % Siegel Z value
        Zm(c) = statsa.zval; % per Matlab
        sig(c) = sig;
        Papprox(c) = pa;
        Pexact(c) = pe;

    else
        % all pairs are identical!
        subs(c) = numel(co);
        nonzero(c) = n;
        nontie(c) = nan;
        Tn(c) = nan;
        Tp(c) = nan;
        Tpr(c) = nan; % rje formula
        Pp(c) = nan;
        Zc(c) = nan; % Conover Z value
        Zs(c) = nan; % Siegel Z value
        Zm(c) = nan; % per Matlab
        sig(c) = nan;
        Papprox(c) = nan;
        Pexact(c) = nan;  
    end

end

%% outputs
output.subs = subs;
output.nonzero = nonzero;
output.nontie = nontie;
output.Tneg = Tn;
output.Tpos = Tp;
output.Tpos_rat = Tpr; % rje formula
output.Ppos = Pp;
output.Z_con = Zc; % Conover Z value
output.Z_sie = Zs; % Siegel Z value
output.Z_mat = Zm; % per Matlab
output.sig = sig;
output.Papprox = Papprox;
output.Pexact = Pexact; 