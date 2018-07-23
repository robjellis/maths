function output = bland_altman_repeated(DATA,SUB,COND,ex_name,co_name)

% method of calculating the bland-altman error bars for repeated measures
% data, as defined by Bland Altman 1999 (sec 5.2)
%
% DATA has two columns: column1 is EXperimental data, column two is COntrol data
%
% COND is any repeated-meausres conditions (e.g., different tempo
% conditions: "uncued", "cued", etc.)
%
% SUB should only be the numbers (1, 2, ... N); gaps are allowed
%
% This program makes use of:
% - anova1, from Matlab stats toolbox
% - plotting function code copied from RJE plot_eda



if size(DATA,1) == size(SUB,1)
    % OK
else
    fprintf(' Error: DATA and SUB are not the same number of rows. \n')
    return
end

if nargin < 3
    COND = [];
end

if nargin < 4

    ex_name = 'EX';
    co_name = 'CO';
end

%% clear NaN manually to be safe
EX = DATA(:,1);
CO = DATA(:,2);

ex_check  = isnan(EX);
co_check  = isnan(CO);
sub_check = isnan(SUB);

all_check = ex_check + co_check + sub_check;
all_check = all_check > 0; % if "1", then it is an NaN

EX  = EX(all_check == 0);
CO  = CO(all_check == 0);
SUB = SUB(all_check == 0);

if numel(COND) == 0
else    
    COND = COND(all_check == 0);
end

sub_id = unique(SUB);
cond_id = unique(COND);

nsub = numel(sub_id);

%% simple correlation
figure(350)
plot(CO,EX,'.')

%% repeated measures stuff

ncond = numel(cond_id);
rm_matrix = nan(nsub,2*ncond); % EX: U, 100, 110; CO: U, 100, 110

for i = 1:nsub
    for j = 1:ncond
        scheck = SUB == sub_id(i);
        ccheck = COND == cond_id(j); % if a cond is missing, the resulting mean will stay at NaN!
        tcheck = scheck + ccheck;
        
        EXk = EX(tcheck == 2);
        COk = CO(tcheck == 2);
        
        rm_matrix(i,j)   = mean(EXk);
        rm_matrix(i,j+3) = mean(COk);
    end
end

%% ANOVAs
displayopt = 'off';

% ANOVA for CO
[p,anovatab,stats] = anova1(CO,SUB,displayopt);
co_mse = cell2mat(anovatab(3,4));

% ANOVA for EX
[p,anovatab,stats] = anova1(EX,SUB,displayopt);
ex_mse = cell2mat(anovatab(3,4));

% how many data points per subject?
pts = stats.n;


%% 1/counts
pts_recip_mn = mean(1 ./ pts);

%% means
mnCO = nan(nsub,1);
mnEX = nan(nsub,1);

sdCO = nan(nsub,1);
sdEX = nan(nsub,1);

for i = 1:nsub
    mnCO(i) = mean(CO(SUB == sub_id(i)));
    mnEX(i) = mean(EX(SUB == sub_id(i)));
    
    sdCO(i) = std(CO(SUB == sub_id(i)));
    sdEX(i) = std(EX(SUB == sub_id(i)));
end

% row mean
mn_row = mean([EX CO],2);

% row diff
diff_row = EX - CO;

% grand mean per subject
mn_grand = mean([mnEX mnCO],2);

% mean differences across subjects
mn_diff = mnEX - mnCO;

% variance of mean differences
var_mn_diff = var(mn_diff);

% mean of mean differences
mn_mn_diff = mean(mn_diff);

%% now the BA plot y-axis variance
varD = var_mn_diff + (1 - pts_recip_mn) * ex_mse + (1 - pts_recip_mn) * co_mse;

% and then we get 1 SD for plotting
stdD = sqrt(varD);

%% and now the error bars
% note: 1.96 is CONSERVATIVE, AND ONLY MAKES SENSE for large sample sizes;
% really, the value should be tinv(.975,df)
% see
% http://en.wikipedia.org/wiki/Inter-rater_reliability#Limits_of_agreement

% so, let's just use "2" instead of 1.96
err_bar = 2 * stdD;
y_plus = mn_mn_diff + err_bar;
y_minus =  mn_mn_diff - err_bar;

% what about the regular error bars for comparison?
err_reg = 2 * std(diff_row);
y_mn_reg = mean(diff_row);
y_plus_reg = y_mn_reg - err_reg;
y_minus_reg = y_mn_reg + err_reg;

%% PLOTS
f199 = figure(199);
set(f199,'Color',[1 1 1])
subplot(1,2,1)
plot(mnCO,sdCO,'.')

subplot(1,2,2)
plot(mnEX,sdEX,'.')

f200 = figure(200);
set(f200,'Color',[1 1 1])

plot(mn_row,diff_row,'Marker','.','MarkerSize',10,'LineStyle','none')
hold on
plot(mn_grand,mn_diff,'MarkerSize',5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0],'Marker','o','LineStyle','none','Color',[1 0 0])

min_x = min(mn_row);
max_x =  max(mn_row);

% plot([min_x max_x],[y_minus y_minus],'g');
% plot([min_x max_x],[mn_mn_diff mn_mn_diff],'g','LineStyle','--');
% plot([min_x max_x],[y_plus y_plus],'g');

% regular error bars
plot([min_x max_x],[y_minus_reg y_minus_reg],'r');
plot([min_x max_x],[y_mn_reg y_mn_reg],'r','LineStyle','--');
plot([min_x max_x],[y_plus_reg y_plus_reg],'r');

title(['Bland-Altman plot' sprintf('\n') '2*SD = ' num2str(err_reg)]);
xlabel(['(' ex_name ' + ' co_name ') / 2'])
ylabel([ ex_name ' - ' co_name ])

hold off

%% Spearman correlations; the p-value is what's important
[rho1 p1] = corr(mnCO,sdCO,'type','Spearman');
[rho2 p2] = corr(mnEX,sdEX,'type','Spearman');
[rho3 p3] = corr(mn_grand,mn_diff,'type','Spearman');

%% OUTPUTS
output.num_sub = nsub;
output.num_pts = size(CO,1);
output.CO_Spearman_p = p1;
output.EX_Spearman_p = p2;
output.BA_Spearman_p = p3;
output.y_mean = mn_mn_diff;
outout.y_std = stdD;
output.y_plus = y_plus;
output.y_minus = y_minus;
output.rm_matrix = rm_matrix;



