function output = ranksum_rje(samp1,samp2,plot_it,do_exact,do_ecdf)

% RJE version of Mann Whitney Wilcoxon 2-sample test
% performs a two-sided rank sum test of the hypothesis
% that two independent samples, in the vectors X and Y, come from
% distributions with equal medians
%
% In particular, this code calcules the "probability of superiority" statistic
%
% calls ranksum.m by Matlab
%
% RJE | 2014.06.23

tic;

if iscell(samp1)
    samp1 = cell2mat(samp1);
end

if iscell(samp2)
    samp2 = cell2mat(samp2);
end

samp1 = samp1(:);
samp2 = samp2(:);

% remove NaN
samp1(isnan(samp1)) = [];
samp2(isnan(samp2)) = [];

if nargin < 3
    plot_it = 0;
end

if nargin < 4
    do_exact = 0;
end

if nargin < 5
    do_ecdf = 0;
end


% do the manual test
n1 = numel(samp1);
n2 = numel(samp2);
nall = n2 + n1;

[ranks, tieadj] = tiedrank([samp1; samp2]);

samp1rank = ranks(1:n1); % sum the ranks of the exerimental group
samp2rank = ranks(n1+1:nall);

% sum the ranks - note: we have already handled ties here since we did
% tiedrank above!
R1 = sum(samp1rank);
R2 = sum(samp2rank); 
%Rtot = nall*(nall+1)/2;

% from wiki:
% http://en.wikipedia.org/wiki/Mann%E2%80%93Whitney_U_test#Calculations
U1 = R1 - n1*(n1+1)/2;
U2 = R2 - n2*(n2+1)/2;
%Utot = n2*n1; % dummy check;

% get the Z score - based directly on the relevant U score
% https://en.wikipedia.org/wiki/Mann%E2%80%93Whitney_U_test#Normal_approximation_and_tie_correction

%Zsamp1_manual = (U1 - n2*n1/2) / sqrt(n1*n2*(n1+n2+1)/12); 
Zsamp2_manual = (U2 - n2*n1/2) / sqrt(n1*n2*(n1+n2+1)/12); % leave this here

% get PS - don't need to do the modification here
%PSsamp1 = U1 / (n1*n2);
PSsamp2 = U2 / (n1*n2);


% add in the Vargha & Delaney (2000) statistic (eq. 14)
%A1 = (R1/n1 - (n1+1)/2)/n2;
A2 = (R2/n2 - (n2+1)/2)/n1;

% compare this with ranksum
[p, ~, stats] = ranksum(samp1,samp2,'method','approximate');
Zsamp1_approx = stats.zval;
Zsamp2_approx = -1 * Zsamp1_approx; % leave this here

if plot_it == 1
    if do_ecdf == 1
        nplot = 3;
    else
        nplot = 2;
    end
    
   figure(90)
   clf
   subplot(1,nplot,1)
   % original values
   maxy = max(max(samp1),max(samp2));
   miny = min(min(samp1),min(samp2));
   range = maxy - miny;
   
   plot(repmat(1,n1,1),samp1,'r.')
   hold on
   plot(repmat(2,n2,1),samp2,'.')
   hold off
   
   axis([0 3 miny-range/100 maxy+range/100])
   ylabel('Original values')
   xlabel('Sample')
   
   subplot(1,nplot,2)
   % plot the ranks
   plot(repmat(1,n1,1),samp1rank,'r.')
   hold on
   plot(repmat(2,n2,1),samp2rank,'.')
   hold off
   axis([0 3 0 nall+1])
   ylabel('Ranks')
   xlabel('Sample')

   if do_ecdf == 1
       % let's also do this in terms of ECDF

       % 1. pool all the data
       data = [samp1(:); samp2(:)];

       % 2. get the ecdf; use RJE program
       f = ecdf_mod(data); % this returns values in same order they were put in

       % 3. pull out the values for sample 1
       f1 = f(1:n1);
       f2 = f(n1+1:end);

       % 4. plot
       subplot(1,nplot,3)
       plot(repmat(1,n1,1),f1,'r.')
       hold on
       plot(repmat(2,n2,1),f2,'.')
       hold off   
       axis([0 3 0 1])
       xlabel('Sample')
       ylabel('ECDF y-axis values')
       

       % 4. take the mean
       f1sum = sum(f1);
       f2sum = sum(f2);
       fsum  = sum(f);

       f1mean = mean(f1);
       f2mean = mean(f2);

       score = 0.5 + (f1mean - f2mean);
   end
   
end

if do_exact == 1
        
    Uobs = min(U2,U1); % this gets samp1mpared against critical value
    
    % for simplicity, let's only worry about n2 == n1
    if n2 ~= n1
        % bad
        error('n2 ~= n1! Uses Z approximation and not samp2act')
    end

    % U critical value (for n2 == n1)
    % these are two-tailed critical values (e.g., Milton 1964, from the "<= 0.25" table
    % 2015.03.16 - RJE samp1nfirmed that if Uobs is <= Ucrit, then we reject the null
    % goes up to N = 30; but Milton only goes up to n = 20 and m = 40
    Ucrit = [nan nan nan nan 2 5 8 13 17 23 30 37 45 55 64 75 87 99 113 127 142 158 175 192 211 230 250 272 294 317];

    if n2 == n1
        H = Uobs <= Ucrit(n2); % <= is key!
    else
        H = 'error';
    end
    % now get the samp2ACT p-value
    pe = ranksum(samp1,samp2,'method','samp2act');
else
    
    pe = 'not evaluated';
end

tocc = toc;

%% outputs
output.numel_samp1 = numel(samp1);
output.numel_unique_samp1 = numel(unique(samp1));
output.numel_samp2 = numel(samp2);
output.numel_unique_samp2 = numel(unique(samp2));
output.num_unique_ranks = numel(unique(ranks));

output.tieadj = tieadj;
output.mean_samp1 = mean(samp1);
output.mean_samp2 = mean(samp2);
%output.R1 = R1;
%output.R2 = R2;
%output.Rtot = Rtot;
output.U1 = U1;
output.U2 = U2;
%output.Utot = Utot;
output.xxx = 'test of samp2 >= samp1';
%output.Zsamp1_manual = Zsamp1_manual;
%output.Zsamp2_manual = Zsamp2_manual;
output.Z_manual = Zsamp2_manual;
%output.Zsamp1_approx = Zsamp1_approx;
%output.Zsamp2_approx = Zsamp2_approx;
output.Z_approx = Zsamp2_approx;
output.Pval_approx = p;
%output.PSsamp1 = PSsamp1;
%output.PSsamp2 = PSsamp2;
output.PS = PSsamp2;
%output.A1 = A1;
%output.A2 = A2;
output.A = A2;

if do_ecdf == 1
   output.f1sum = f1sum;
   output.f2sum = f2sum;
   output.f1mn = f1mean;
   output.f2mn = f2mean;
   output.f_score = score;
end
output.duration = tocc;


if do_exact == 1
    output.Hsamp2act   = H; % 0 = not sig.; 1 = sig.
    output.Psamp2act = pe;
end
