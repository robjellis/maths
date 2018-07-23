function out = propsup_boot(x1,x2,expN,iter,plot_PS,do_MW)

tic;
% manually derive the "proportion of superiority" statistic

if nargin < 3
    expN = 20;
    iter = 5000;
end

if nargin < 5
    plot_PS = 0;
    do_MW = 0;
end

% 1. How many?
n1 = numel(x1);
n2 = numel(x2);

% 2. get indices with replacement (i.e., bootstrap) from samp1 and samp2
i1 = randint(expN,iter,[1 n1]);
i2 = randint(expN,iter,[1 n2]);

% 3. match up with values
y1 = x1(i1);
y2 = x2(i2);

% 5. Count each COLUMN separately

% 5a. count superior
sup = sum(y1 > y2);

% 5b. count ties
ties = sum(y1 == y2) / 2; % count ties at 50%

% 5c. get the PS
PS = (sup + ties) / expN;

% 6. get statistics
PS5 = prctile(PS,5);
PS50 = prctile(PS,50);
PSMN = mean(PS);

% how many are >= 2/3
hi = sum(PS >= 2/3);

lo = sum(PS <= 1/3);

prop_hi = hi / iter;
prop_lo = lo / iter;

EC = (hi + -1*lo) / iter; % "experiment consistency"

% do MW PS calc?
if do_MW == 1
    rs = ranksum_rje(x1,x2);
    MW_PS = rs.PSsamp1;
else
    MW_PS = NaN;
end

tocc = toc;

% 7. Now get distribution of PS
if plot_PS == 1
    figure(100)
    ecdf(PS)
    
end

%% output
out.samp_size = expN;
out.iter = iter;
out.PS5 = PS5;
out.PS50 = PS50;
out.PSMN = PSMN;
out.MW_PS = MW_PS;
out.EC    = EC;
out.prop_hi = prop_hi;
out.prop_lo = prop_lo;
out.duration = tocc;

