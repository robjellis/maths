function out = ranksum_1samp(varargin)

% RJE modification of MW test to work with a *single* sample of data
% note: currently formatted to test against normal distribution, but other
% distribution shapes could easily be added
%
% assumes data in samp1 is a single vector
%
% Jan 12 2018

samp1 = varargin{1};
samp1 = samp1(:);

N = numel(samp1);

if nargin == 1
    % determine bounds for simulated
    lv = (1 / N) / 2;
    uv = 1 - lv;

    samp2 = norminv(linspace(lv,uv,N)); % can't use 0 and 1 as bounds since those yield -Inf and +Inf, respectively
	samp2 = samp2(:);
else
    samp2 = varargin{2};
    samp2 = samp2(:);
end

% do the manual test
n1 = numel(samp1);
n2 = numel(samp2);

[ranks, tieadj] = tiedrank([samp1; samp2]);

% how many tied values?
num_ties = numel(ranks) - numel(unique(ranks));

samp1rank = ranks(1:n1); % sum the ranks of the exerimental group

% sum the ranks - note: we have already handled ties here since we did tiedrank above!
R1 = sum(samp1rank);

% from wiki:
% http://en.wikipedia.org/wiki/Mann%E2%80%93Whitney_U_test#Calculations
U1 = R1 - n1*(n1+1)/2;

% get the Z score (normal approximation) - based directly on the relevant U score
% https://en.wikipedia.org/wiki/Mann%E2%80%93Whitney_U_test#Normal_approximation_and_tie_correction
Zman = (U1 - n2*n1/2) / sqrt(n1*n2*(n1+n2+1)/12); 

% 2-tailed p-value
if Zman < 0
    Pman = normcdf(Zman) * 2;    
else
    Pman = normcdf(-Zman) * 2; % we get more precise p-values when estimating the low end of the CDF
end


% get PS - don't need to do the modification here
PSsamp1 = U1 / (n1*n2);

out.Zman     = Zman;
out.Pman     = Pman;
out.PS1      = PSsamp1;
out.num_ties = num_ties;
out.samp2    = samp2;