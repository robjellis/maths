function output = ba_repeated(STRUC,DATA,outlier_meth,outlier_thr,ex_name,co_name,ba_mod,do_sim)

%
% output = ba_repeated(STRUC,DATA,outlier_thr,ex_name,co_name,ba_mod)
% 
% method of calculating the bland-altman error bars for repeated measures
% data, as defined by Bland Altman 1999 (sec 5.2)
%
%
% STRUC is an [N x 4] matrix with columns in this order:
%  1. INCL  is the single means of removing any cases; only INCL == 1 will be retained
%  2. GROUP is to distinguish between patients (GROUP = 1) and controls (GROUP = 2) (use a single [] to ignore)
%  3. SUB   should only integers (1, 2, ...); gaps are allowed
%  4. COND  is any repeated-meausres conditions (e.g., different tempo conditions: "uncued", "cued", etc.)
%
% DATA is comprised of two columns: 
%  1. EX is experimental (novel) measure 
%  2. CO is control (standard) measure

% outlier_meth is:
%    'n'   = none; just plot y-axis mean and y-axis 2*SD separately for each group
%    'sd'  = outlier detection mode: outliers defined as falling outside
%            interval [mean - m * SD, mean + m * SD], where m is "outlier_thr"
%    'mad' = robust outlier detection: outliers defined as falling outside
%            interval [median - m * MAD, median + m * MAD], where m is "outlier_thr"
%
% This program makes use of:
% - anova1, from Matlab stats toolbox
% - plotting function code copied from RJE plot_eda

EX = DATA(:,1);
CO = DATA(:,2);


%% color or black and white?
clr = 0;


if isempty(STRUC)
    INCL  = ones(size(CO)); % retain all cases
    GROUP = ones(size(CO)); % just one group
    SUB   = 1:size(CO,1);
    COND  = ones(size(CO));
else
    INCL  = STRUC(:,1);
    GROUP = STRUC(:,2);
    SUB   = STRUC(:,3);
    COND  = STRUC(:,4);
end

if nargin < 3
    outlier_meth = 'n'; % no outlier removal
    outlier_thr = [];
end

if nargin < 5
    ex_name = 'EX';
    co_name = 'CO';
end

if nargin < 7
    ba_mod = 0; % if "1", then will just use CO values on X-axis (and not average of EX and CO)
    do_sim = 0;
end


%% clear NaNs based on INCL

EX    = EX(INCL == 1);
CO    = CO(INCL == 1);
SUB   = SUB(INCL == 1);
COND  = COND(INCL == 1);
GROUP = GROUP(INCL == 1);

% further reduce in case we only have one condition
EX    = EX(COND > 0);
CO    = CO(COND > 0);
SUB    = SUB(COND > 0);
GROUP   = GROUP(COND > 0);
COND    = COND(COND > 0);

sub_id = unique(SUB);
cond_id = unique(COND);
group_id = unique(GROUP);

nsub = numel(sub_id);
ncond = numel(cond_id);
ngroup = numel(group_id);


%% output variables
% this only is useful when ncond == 1
if ncond == 1
    y_mean_all = nan(1,ncond);
    y_std_all  = nan(1,ncond);
    CO_mean_all = nan(1,ncond);
    CO_std_all  = nan(1,ncond);
    
end

%% simple correlation
figure(350)
plot(CO,EX,'.')

% the original data
EXo = EX;
COo = CO;
GROUPo = GROUP;

% a copy to save prior to the group loop
EXf   = EX;
COf   = CO;
SUBf  = SUB;
CONDf = COND;
GROUPf = GROUP;


%% IDENTIFY OUTLIERS ACROSS THE FULL DATASET
% This is an optional step
% Do this pooled over all subjects and all conditions to get the most
% stable estimate possible


if strcmp(outlier_meth,'n')
    % no outlier removal
    EX_out = [];
    CO_out = [];
    SUB_out = [];
elseif strcmp(outlier_meth,'mad') || strcmp(outlier_meth,'sd')
    dev_err = EX - CO; % all values pooled together over groups and conditions
    out_plot = 1;
    outID = outlier_rje(dev_err,outlier_meth,outlier_thr,out_plot);
    
    % now actually get the outlier points
    % these are pooled over groups but it doesn't matter; we plot them all the same color
    EX_out  = EXf(outID.inds);
    CO_out  = COf(outID.inds);
    SUB_out = SUBf(outID.inds);
    
    % clean it up
    EXf(outID.inds)   = [];
    COf(outID.inds)   = [];
    SUBf(outID.inds)  = [];
    CONDf(outID.inds) = [];
    GROUPf(outID.inds) = [];

end

cond_std_max = nan(ngroup,1);

%% group-level loop    
for g = 1:ngroup

    % now redefine the variables
    EX    = EXf(GROUPf == group_id(g));
    CO    = COf(GROUPf == group_id(g));
    SUB   = SUBf(GROUPf == group_id(g));
    COND  = CONDf(GROUPf == group_id(g));
    
    % a copy for plotting purposes
    EXp = EXo(GROUPo == group_id(g));
    COp = COo(GROUPo == group_id(g));

    ba_y = EXp - COp;
    ba_y_out = EX_out - CO_out;

    if ba_mod == 1
        ba_x = COp;
        
        ba_x_out = CO_out;
    else
        ba_x = (EXp + COp) / 2;
        ba_x_out = (EX_out + CO_out) / 2;
    end
    
    % it may be the case that certain *data* conditions have NaN; we need to
    % remove these cases so it doesn't cause problems for ANOVA
    
    exN = isnan(EX);
    coN = isnan(CO);
    
    allN = max(exN,coN);
    
    EX   = EX(allN == 0);
    CO   = CO(allN == 0);
    SUB  = SUB(allN == 0);
    COND = COND(allN == 0);
    
    % set the colors
    if clr == 1
        if g == 1
            col = 'r';
            width = 2.0;
        elseif g == 2
            col = 'b';
            width = 1.0;
        elseif g == 99
            col = 'm'; % merged
            width = 1.75;
        end
        col_out = 'm'; % for marking outliers
        
    elseif clr == 0
        if g == 1
            col = 'k'; % black
            width = 2.0;
        elseif g == 2
            col = [.5 .5 .5]; % gray
            width = 1.0;
        elseif g == 99
            col = 'm'; % merged
            width = 1.75;
        end       
        col_out = [.25 .25 .25];
    end
    
    
    %% repeated measures stuff
    
    if ncond > 1
        rm_matrix = nan(nsub,1 + 2*ncond); % EX: U, 100, 110; CO: U, 100, 110
        for i = 1:nsub

            scheck = SUB == sub_id(i); % only do this for the subject IDs we actually use to save time
            for j = 1:ncond         
                ccheck = COND == cond_id(j); % if a cond is missing, the resulting mean will stay at NaN!
                tcheck = scheck + ccheck;

                EXk = EX(tcheck == 2);
                COk = CO(tcheck == 2);

                rm_matrix(sub_id(i),1) = sub_id(i);
                rm_matrix(sub_id(i),1 + cond_id(j))       = mean(EXk);
                rm_matrix(sub_id(i),1 + cond_id(j)+ncond) = mean(COk);
            end
        end
        
        % RJE: get the 2*STD for each condition and then take the max; do
        % this separately for each group
        
        mat_ex = rm_matrix(:,2:ncond+1);
        mat_co = rm_matrix(:,ncond+2:end);
        rm_matrix_diff = mat_ex - mat_co;   % EX minus control
        
        rm_matrix_diff_std = nanstd(rm_matrix_diff);
        
        cond_std_max(g) = max(rm_matrix_diff_std);
        
        figure(40)
        if ncond == 1
            clf
        end
        if g == 1
            clf
        else
            hold on
        end

        subplot(2,1,1)
        dev_err = EX-CO;
        plot(SUB,dev_err,'MarkerSize',10,'MarkerFaceColor',col,'MarkerEdgeColor',col,'Marker','.','LineStyle','none','Color',col)
        hold on

        % identify the outliers
        if isempty(outlier_thr) == 0
            if g == ngroup

            % we already defined the parameters for outlier detection above
            %plot(SUB_out,EX_out-CO_out,'mo')
            end
        end

        if g == ngroup
            hold off
        end
        ylabel([ ex_name ' - ' co_name ])

        subplot(2,1,2)

        % within-subject standard deviation
        std_flag = 0; % means use N - 1
        sub_std = nanstd((rm_matrix(:,2:ncond+1) - rm_matrix(:,ncond+2:end)) ,std_flag,2);

        plot(rm_matrix(:,1),sub_std,'MarkerSize',10,'MarkerFaceColor',col,'MarkerEdgeColor',col,'Marker','.','LineStyle','none','Color',col)
        if g < ngroup
            hold on
        else
            hold off
        end
        ylabel(['std(' ex_name ' - ' co_name ')'])



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
            % use NaN mean in case we have a condition with NaN
            mnCO(i) = nanmean(CO(SUB == sub_id(i)));
            mnEX(i) = nanmean(EX(SUB == sub_id(i)));

            sdCO(i) = nanstd(CO(SUB == sub_id(i)));
            sdEX(i) = nanstd(EX(SUB == sub_id(i)));
        end

        % BA individual x-values
    %     if ba_mod == 1
    %         mn_row = CO; % just use ground truth data
    %     else
    %         mn_row = nanmean([EX CO],2);
    %     end

        % BA individual y-values
        diff_row = EX - CO;

        % grand mean per subject
        if ba_mod == 1
            mn_grand = mnCO;
        else
            mn_grand = mean([mnEX mnCO],2);
        end

        % mean differences across subjects
        mn_diff = mnEX - mnCO;

        % variance of mean differences
        var_mn_diff = nanvar(mn_diff);

        % mean of mean differences
        y_mean = nanmean(mn_diff);

        %% now the BA plot y-axis variance
        varD = var_mn_diff + (1 - pts_recip_mn) * ex_mse + (1 - pts_recip_mn) * co_mse;

        % and then we get 1 SD for plotting
        stdD = sqrt(varD);
        y_std = stdD; % just for simplicity

        %% and now the error bars
        % note: 1.96 is CONSERVATIVE, AND ONLY MAKES SENSE for large sample sizes;
        % really, the value should be tinv(.975,df)
        % see
        % http://en.wikipedia.org/wiki/Inter-rater_reliability#Limits_of_agreement

        % so, let's just use "2" instead of 1.96
        err_bar = 2 * stdD;
        err_bar_all(g) = err_bar;
        y_plus = y_mean + err_bar;
        y_minus =  y_mean - err_bar;

        % what about the regular error bars for comparison?
        % 2014.10.14 - we now do this as just taking the largest error across any condition

        err_reg = 2 * cond_std_max(g);
        err_reg_all(g) = err_reg;
        y_mn_reg = mean(diff_row);
        y_plus_reg = y_mn_reg - err_reg;
        y_minus_reg = y_mn_reg + err_reg;
    elseif ncond == 1
        % only one cond
        y_mean = nanmean(EX-CO);
        y_std = nanstd(EX-CO);
        err_bar = 2 * y_std;
        y_plus = y_mean + err_bar;
        y_minus = y_mean - err_bar;    
        
        
        CO_mean  = nanmean(ba_x);
        CO_std   = nanstd(ba_x);
        
        % we force this to just be the CO data
        CO_mean_all(g) = nanmean(CO);
        CO_std_all(g) = nanstd(CO);
        
        y_mean_all(g) = y_mean;
        y_std_all(g) = y_std;
    end

    %% PLOTS
%     f1 = figure(199);
%     set(f1,'Color',[1 1 1])
%     subplot(1,2,1)
%     plot(mnCO,sdCO,'.')
% 
%     subplot(1,2,2)
%     plot(mnEX,sdEX,'.')

%     if g == 1
%         f2 = figure(201);
%     elseif g == 2
%         f2 = figure(202);
%     end
    
    if g == 1
        f2 = figure(200);
        clf
        pause(.2)
        set(f2,'Color',[1 1 1])
    else
        figure(200)
        hold on % so we plot the ohter data on top
    end
    

    plot(ba_x,ba_y, 'MarkerSize',2, 'MarkerFaceColor',col,'MarkerEdgeColor',col,'Marker','o','LineStyle','none','Color',col)
    hold on
    
    if ncond > 1
        plot(mn_grand,mn_diff,'MarkerSize',5, 'MarkerFaceColor',col,'MarkerEdgeColor',col,'Marker','o','LineStyle','none','Color',col)
        min_x = min(mn_grand);
        max_x =  max(mn_grand);
        
        plot_reg = 0;
        if plot_reg == 1
            plot([min_x max_x],[y_minus y_minus],'g');
            plot([min_x max_x],[y_mean y_mean],'g','LineStyle','--');
            plot([min_x max_x],[y_plus y_plus],'g');
        end

        % regular error bars
        if strcmp(outlier_meth,'n')
            plot([min_x max_x],[y_minus_reg y_minus_reg],'Color',col,'LineWidth',width);
            plot([min_x max_x],[y_mn_reg y_mn_reg],      'Color',col,'LineWidth',width,'LineStyle','--');
            plot([min_x max_x],[y_plus_reg y_plus_reg],  'Color',col,'LineWidth',width);
        end

    else
        min_x = min(ba_x);
        max_x = max(ba_x);

        
        % do the standard 95% CI, using Cumming and Fitch 2005 as a
        % reference; for a single group of data,   95% CI = tinv(.975,df) * STD / sqrt(n)
        % actually, that will be misleading because we are using average of
        % EX and CO on the x-axis!
        
        % how many non-NaN values are there?
        num_val = sum(isnan(ba_x)==0);
        

        err_width = 1;
        
        if err_width == 1
            % just the range of the non-outliers
            if ba_mod == 1
                x_plus = max(CO);
                x_minus = min(CO);
            elseif ba_mod == 0
                x_plus = max((CO+EX)/2);
                x_minus = min((CO+EX)/2);
            end
        elseif err_width == 2
            % 2 * STD
            x_plus  = CO_mean + 2 * CO_std;
            x_minus = CO_mean - 2 * CO_std; 
        elseif err_width == 3
            % CI
            x_plus  = CO_mean + 2 * CO_std / sqrt(num_val);
            x_minus = CO_mean - 2 * CO_std / sqrt(num_val);
        
        end
        
        if strcmp(outlier_meth,'n')
            plot([x_minus x_plus],[y_minus y_minus],'Color',col,'LineWidth',width);
            plot([x_minus x_plus],[y_mean y_mean],  'Color',col,'LineWidth',width,'LineStyle','--');
            plot([x_minus x_plus],[y_plus y_plus],  'Color',col,'LineWidth',width);
        elseif strcmp(outlier_meth,'mad') || strcmp(outlier_meth,'sd')
            % we only do this once
            if g == ngroup
                xmin_full = min((EXo+COo)/2);
                xmax_full = max((EXo+COo)/2);
                plot([xmin_full xmax_full],[outID.thr_lo outID.thr_lo],'k','LineStyle','--');
                plot([xmin_full xmax_full],[outID.thr_hi outID.thr_hi],'k','LineStyle','--');  
            end
        end
        
        % add the x-axis mean
        if ba_mod == 1
            plot(CO_mean,y_minus,'MarkerSize',5, 'MarkerFaceColor',col,'MarkerEdgeColor',col,'Marker','square','LineStyle','none','Color',col)
            plot(CO_mean,y_plus ,'MarkerSize',5, 'MarkerFaceColor',col,'MarkerEdgeColor',col,'Marker','square','LineStyle','none','Color',col)
        end
        
    end
    
    % mark the outliers
    plot(ba_x_out,ba_y_out,'MarkerSize',8,'Marker','o','LineStyle','none','Color',col_out)
    
%     % adjust th axis scale
%     xrange = max(ba_x) - min(ba_x);
%     yrange = max(ba_y) - min(ba_y);
%     
%     div = 4;
%     axis([min_x-xrange/div    max_x+xrange/div    min(ba_y)-yrange/div   max(ba_y) + yrange/div])
    

    if ngroup == 1
        title(['Bland-Altman plot' sprintf('\n') '2*SD = ' num2str(err_bar)]);
    else
        %title('Bland-Altman plot');
    end

    if ba_mod == 1
        xlabel(co_name)
    else
        xlabel(['(' ex_name ' + ' co_name ') / 2'])
    end

    ylabel([ ex_name ' - ' co_name ])


    hold off


end % end of group loop

%% Spearman correlations; the p-value is what's important

if ncond > 1
[rho1 p1] = corr(mnCO,sdCO,'type','Spearman');
[rho2 p2] = corr(mnEX,sdEX,'type','Spearman');
[rho3 p3] = corr(mn_grand,mn_diff,'type','Spearman');
else
p1 = [];
p2 = [];
[rho3 p3] = corr(CO,EX,'type','Spearman');    
end


%% simulation
% how often do we get a significant result between Group1 and Group2 when
% we use the CO device versus the EX device?
if do_sim == 1 && ngroup == 2
   fprintf('\n Performing simulation ... \n')
   
   n = 15;
   iter = 1000;
   
   MWsig = nan(iter,2); % to store significance (1 = significant); col1 is for CO, col2 is for EX
   
   for i = 1:iter
       % it is possible to get negative values, so we create 2N values and
       % then just take the first N; this adds a bit more variability which is OK
       
       g1 = randn(2*n,1) * CO_std_all(1) + CO_mean_all(1);
       g2 = randn(2*n,1) * CO_std_all(2) + CO_mean_all(2);

       
       % now, importantly, we ADD IN the error induced by the device  
       g1co = g1;
       g1ex = g1 + randn(2*n,1) * y_std_all(1) + y_mean_all(1);
       
       g2co = g2;
       g2ex = g2 + randn(2*n,1) * y_std_all(2) + y_mean_all(2);
       
       % restrict to positive values
       g1neg = max([g1co<0 g1ex<0],[],2);
       g2neg = max([g2co<0 g2ex<0],[],2);
       
       % eliminate negative values? only if original data values are all
       % positive
       
       if min(CO) > 0
           g1co(g1neg==1) = [];
           g1ex(g1neg==1) = [];

           g2co(g2neg==1) = [];
           g2ex(g2neg==1) = [];
       end
       
       
       % take the first N
       g1co = g1co(1:n);
       g1ex = g1ex(1:n);
       
       g2co = g2co(1:n);
       g2ex = g2ex(1:n);

        % do the non-parametric Mann Whitney test for robust significance
        % http://www.mathworks.com/help/stats/ranksum.html
        % g1 and g2 can have different numbers of values

        [p hco] = ranksum(g1co,g2co); 
        [p hex] = ranksum(g1ex,g2ex); 
        
        MWsig(i,1) = hco;
        MWsig(i,2) = hex;
   end
   
   
   MWsig = 100 * sum(MWsig) / iter;
end
%% OUTPUTS
output.num_sub = nsub;
output.num_pts = size(CO,1);
output.CO_Spearman_p = p1;
output.EX_Spearman_p = p2;
output.BA_Spearman_p = p3;



if ncond > 1
    output.y_err_bar_BA = err_bar_all;
    output.y_err_bar_rje = err_reg_all;
    output.rm_matrix = rm_matrix;
elseif ncond == 1
    output.y_mean = y_mean_all;
    output.y_2std = 2 * y_std_all;
    output.y_plus = y_plus;
    output.y_minus = y_minus;

    output.CO_mean = CO_mean_all;
    output.CO_std = CO_std_all;
end

if do_sim == 1
   output.sim_iter = iter;
   output.sig_CO = MWsig(1);
   output.sig_EX = MWsig(2);
end


