function output = lr_calc(group_code,data_all,kde_meth,plot_it)

% output = lr_calc(group_code,data_all,kde_meth,plot_it)
%
% calulates ikelihood ratios using the definitions by Choi et al., 1988:
% specifically, LR(x,y)
%
% group_code: an [N x 1] of 1s (= healthy controls, no disease: "HC") or 2s (= patients: "PT"); 0 means ignore the case
% data: a column of observed statistic values
%
% data_all: an [N x M] matrix, with each M being a distinct variables / outcome measure
%
% Choi (1988). Slopes of a receiver operating characteristic curve and
% likelihood ratios for a diagnostic test. American journal of
% epidemiology, vol 148(11), p. 1127-1132.
%
% rje | 2013.09.18

if nargin < 3
    plot_it = 0;
end

% check size
if size(group_code,1) == size(data_all,1)
    % OK
else
    fprintf('\n Warning: the number of subject cases in group_code does not match the number of subject cases in data_all.\n\n')
    return
end

% get rid of NaNs
keep_inds = isnan(group_code) == 0; % 1 means keep
 
data_all   = data_all(keep_inds,:);
group_code = group_code(keep_inds);

% how many HC and PT?

num_HC = sum(group_code == 1);
num_PT = sum(group_code == 2);

% report it
fprintf(['\n Read in ' num2str(num_HC) ' controls and ' num2str(num_PT) ' patients. \n'])
% how to do the LR

LRlook = 0.5:0.05:0.9;

% set up output variables
ncols = size(data_all,2);

LR_xval = nan(1,ncols);
LR_look = LRlook';
LR_all = nan(numel(LRlook),ncols);
LR_max = nan(1,ncols);
cdf_vertdis_xval = nan(1,ncols);
cdf_vertdis_yval = nan(1,ncols);
cdf_propcor_xmin = nan(1,ncols);
cdf_propcor_ymin = nan(1,ncols);
cdf_propcor_zval = nan(1,ncols);

%% loop

for m = 1:ncols
    
    data = data_all(:,m);
    
    % separate the variables
    data_hc = data(group_code == 1);
    data_pt = data(group_code == 2);
    
    % number of subjects

    %% KS test
    [h p k2stat] = kstest2(data_hc,data_pt,.05,'smaller');


    %% KDE estimation
    % http://www.mathworks.com/matlabcentral/fileexchange/14034-kernel-density-estimator

    nmesh = 2^12; % default = 2^12

    xinds = 1:nmesh;
    
    % we need a COMMON scale for the data
    range = max(data) - min(data);
    MIN = min(data) - range/4; 
    MAX = max(data) + range/4;

    if strcmp(kde_meth,'m')
        % matlab
        xvals = linspace(MIN,MAX,nmesh);
        X_hc = xvals;
        X_pt = xvals;

        [D_hc]   = ksdensity(data_hc,X_hc,'function','pdf');
        [CDF_hc] = ksdensity(data_hc,X_hc,'function','cdf');

        [D_pt]   = ksdensity(data_pt,X_pt,'function','pdf');
        [CDF_pt] = ksdensity(data_pt,X_pt,'function','cdf');

    elseif strcmp(kde_meth,'b')
        % botev
        [bw_hc,D_hc,X_hc,CDF_hc] = kde(data_hc,nmesh,MIN,MAX); % control
        [bw_pt,D_pt,X_pt,CDF_pt] = kde(data_pt,nmesh,MIN,MAX); % patint
    end

    % make sure we are bounded correctly
    D_hc(D_hc < 0) = 0;
    D_pt(D_pt < 0) = 0;
    
    CDF_hc(CDF_hc < 0) = 0;
    CDF_hc(CDF_hc > 1) = 1;

    CDF_pt(CDF_pt < 0) = 0;
    CDF_pt(CDF_pt > 1) = 1;

    %% RJE calculations

    % % PDF: HC> PT
    % pdf_HC_min_PT = D_hc - D_pt;
    % 
    %     % what is the first zero crossing?
    %     pdf_xval_cross = min(X_pt(pdf_HC_min_PT <= 0));
    % 
    %     figure(301)
    %     plot(X_pt,pdf_HC_min_PT,'b')
    %     hold on
    %     plot([MIN MAX], [0 0],'r')
    %     hold off

    % PDF: max(HC,PT)
    pdf_min = min(D_hc,D_pt); % at every point

    % CDF: vertical difference between HC and PT; similar to KS test
    HC_min_PT = CDF_hc - CDF_pt; % sign matters: *only* positive sign is useful in a diagnostic sense

    cdf_ydis_max = max(HC_min_PT);
    cdf_xval_max = min(X_pt(HC_min_PT == max(HC_min_PT)));

    % x-value that maximizes both true positives and true negatives
        % this version may not be necessary
        cdf_propcor_min = min((1 - CDF_pt),CDF_hc); % 1 - CDF is the complementary CDF (CCDF)

        % this is better: weighted version of true pos. and true neg.
        pos_w = 1;
        neg_w = 1;
        cdf_propcor_weight = (pos_w*(1 - CDF_pt) + neg_w*CDF_hc) / (pos_w + neg_w);
        
        
        % RJE: there may be values slightly less than 1.0 (due to aspects of KDE estimation), so we need to adjust
        if max(cdf_propcor_min) == 1
           cdf_propcor_min(cdf_propcor_min >= .9999) = 1; % leave this rule here unless we figure out the cause
        end

        % ind_max = cdf_propcor_min == max(cdf_propcor_min)
        % zz = [X_pt(:) cdf_propcor_min(:) ind_max(:)]

        cdf_propcor_ymin = max(cdf_propcor_min);
        cdf_propcor_xmin = X_pt(cdf_propcor_min == cdf_propcor_ymin);

        % there may be multiple values, so just take the median
        cdf_propcor_xmin = median(cdf_propcor_xmin);

        % now to get this to a "useful" scale, take the 
        % Normal inverse cumulative distribution function (i.e., find the corresponding Z-value)

        % note: maybe this is not useful, because the "raw" values are at least
        % interpretable; Z-values are not!

        cdf_propcor_zval = norminv(max(cdf_propcor_ymin));

        % of course, we have to correct for P = 1.0 which means Z = Inf

        if cdf_propcor_zval == Inf;
           cdf_propcor_zval = 10; % still way higher than any other possible Z-value, which is ~ 8.21
        end


    %% LR calculation
    % LR(x,Inf); equivalent to saying "likelihood of disease given a score of x or higher"

    % note: we constrain the "look" range to x-values between (1) CDF_pt <= .01 and (2) CDF_hc >= .99

    thr = .05;

    x_low_pt  = max(X_pt(CDF_pt <= thr));
    x_high_hc = min(X_hc(CDF_hc >= (1-thr)));


    % are we OK?
    if x_low_pt > x_high_hc
        % the PT group is much higher than HC; the test has very high discrimination!
        LR_xval = (x_low_pt + x_high_hc) / 2;
        LR_max = 1000; % just to return a value

    else

        % now get the index values to look at
        x1 = find(X_pt >= x_low_pt);
        x2 = find(X_hc <= x_high_hc);

        x_look_inds = intersect(x1,x2);
        x_look_vals = X_pt(x_look_inds);
        
        % 2013-09-19: update - since we are more interested in accurately
        % detecting PT group, let's just look at LR for various values on
        % the CDF for PT - e.g., 50:10:90; 
        
        for y = 1:numel(LRlook)
           xind(y) = max(xinds((1-CDF_pt)>LRlook(y))); % these will be index positions
        end

        xlook = X_pt(xind);
        
        % this is simply done in terms of the CDF
        LR = (1 - CDF_pt(xind)) ./ (1 - CDF_hc(xind));

        % find the max LR
        LR_xval  = min(x_look_vals(LR == max(LR)));
        LR_max  = max(LR);

    end


    %% plots

    if plot_it == 1
       figure(500)
       subplot(2,7,1)
       plot(group_code,data,'.')
       xlabel('Condition')
       ylabel('Statistic')
       axis([0 3 MIN MAX])

       % PDF
       subplot(2,7,[2 3])
       plot(X_hc,D_hc,'g')
       hold on
       plot(X_pt,D_pt,'r')
       hold off
       xlabel('Statistic')
       ylabel('Estimated probability')
       title('PDF')
       ymax = max(max(D_hc),max(D_pt));
       axis([MIN MAX 0 ymax*1.1])

       subplot(2,7,[9 10])
       plot(X_hc,pdf_min,'b')
       xlabel('Statistic')
       ylabel('min(HC,PT)')
       axis([MIN MAX 0 max(pdf_min)*1.1])

       % CCDF
       subplot(2,7,[4 5])
       plot(X_hc,1 - CDF_hc,'g')
       hold on
       plot(X_pt,1 - CDF_pt,'r')
       %plot(X_pt,HC_min_PT,'b')
       hold off
       xlabel('Statistic')
       ylabel('1 - Pr(D<=X)')
       title('CCDF')
       axis([MIN MAX 0 1])   

       subplot(2,7,[11 12])
       plot(X_hc,CDF_hc,'g')
       hold on
       plot(X_pt,1-CDF_pt,'r') % CCDF
       %plot(X_hc,cdf_propcor_min,'b')
       plot(X_hc,cdf_propcor_weight,'b')
       hold off
       xlabel('Statistic')
       ylabel('Maximum true neg. (HC) and true pos. (PT)')   
       axis([MIN MAX 0 1])

       % only plot LR if it is valid
       if x_low_pt < x_high_hc
           % LR
           subplot(2,7,[6 7])
           plot(0,0) % just to clear it
           %plot(x_look_vals,LR,'b')
           plot(xlook,LR,'b')
           xlabel('Statistic threshold')
           ylabel('LR(x,Inf)')
           axis([MIN MAX 0 max(LR)*1.1]) 
       end

       if ncols > 1
           pause(0.5)
       end
    end

    %% outputs
    output.num_HC = num_HC;
    output.num_PT = num_PT;
    output.mesh_points = nmesh;
    output.LR_xval(m) = LR_xval;
    output.LR_look = LR_look;
    output.LR_all(:,m) = LR(:);
    output.LR_max(m) = LR_max;
    output.cdf_vertdis_xval(m) = cdf_xval_max;
    output.cdf_vertdis_yval(m) = cdf_ydis_max;
    output.cdf_propcor_xmin(m) = cdf_propcor_xmin;
    output.cdf_propcor_ymin(m) = cdf_propcor_ymin;
    output.cdf_propcor_zval(m) = cdf_propcor_zval;

end


