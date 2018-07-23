function output = lr_calc(group_code,data,nprc,plot_it)

% output = lr_calc(group_code,data,nprc,plot_it)
%
% calulates positive likelihood ratios (LR+) and returns the threshold (for the data) and
% the maximum observed LR+
%
% group_code: a column of 1s (= control) or 2s (= patient); 0 means ignore the case
% data: a column of observed statistic values
%
% rje | 2013.09.18

if nargin < 3
    nprc = 50;
end

if nargin < 4
    plot_it = 0;
end

if sum(size(group_code) == size(data)) == 2
    % OK
    group_code = group_code(:);
    data = data(:);
else
    return
end

% get rid of zeros
data(group_code == 0) = [];
group_code(group_code == 0) = [];

% separate the variables
data_g1 = data(group_code == 1);
data_g2 = data(group_code == 2);

% create the percentiles
minp = 10;
maxp = 100 - minp;

prc_vals = linspace(minp,maxp,nprc);

% get percentiles from each condition *separately*
prc_thr_d1 = prctile_nist(data_g1,prc_vals);
prc_thr_d2 = prctile_nist(data_g2,prc_vals);

% merge them
prc_thr_all = [prc_thr_d1; prc_thr_d2];
prc_thr_all = sort(unique(prc_thr_all));

nprc_all = numel(prc_thr_all);

% set up holding variables
sens_all = nan(nprc_all,1);
spec_all = nan(nprc_all,1);



%% loop

for p = 1:nprc_all
   
    % sensitivity
    sens_all(p) = sum(data_g2 > prc_thr_all(p)) / numel(data_g2); % 0 to 1 value
    
    % specificity
    spec_all(p) = sum(data_g1 < prc_thr_all(p)) / numel(data_g1); % 0 to 1 value
    
end


%% find the max LR+

lr_plus = sens_all ./ (1 - spec_all);


%% plots

if plot_it == 1
   figure(300)
   subplot(1,4,1)
   plot(group_code,data,'.')
   hold on
   plot(zeros(numel(prc_thr_all,1)),prc_thr_all,'r+')
   xlabel('Condition')
   ylabel('Statistic')
   
   range = max(data) - min(data);
   minn = min(data) - range/10;
   maxx = max(data) + range/10;
   axis([0 3 minn maxx])
   
   subplot(1,4,[2 3])
   plot(prc_thr_all,sens_all,'r')
   hold on
   plot(prc_thr_all,spec_all,'b')
   hold off
   xlabel('Statistic threshold')
   ylabel('Proportion')
   axis([minn maxx 0 1])
   
   subplot(1,4,4)
   plot(prc_thr_all,lr_plus,'g')
   xlabel('Statistic threshold')
   ylabel('LR+')
   axis([minn maxx 0 1.1 * max(lr_plus)]) 
end

%% outputs

output.lr_max = lr_plus(lr_plus == max(lr_plus));
output.thr_max = prc_thr_all(lr_plus == max(lr_plus));


