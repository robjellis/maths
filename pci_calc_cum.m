function [res] = pci_calc_cum(Xmat,plot_phi)

%
% PCI calculation by Plotnik et al. 2007, with additional calculations by RJE
%
% allows *multiple* runs within Xmat to be identified; phi calculated; and
% then PCI calculated. We do this "improved" method for both the Hausdorff
% PCI method, and RJE merged method, just to level the playing field
% 
% Note: PCI is only useful for gait data, as it gives information about
% alternation of a two-element process (L vs R). For tempo data or HR data,
% use the more general td_calcs.m by RJE.
%
% Xmat is a four-column matrix that is output from td_outlier_detect by
% RJE. col1 = timestamps; col2 = channel (L = 1, R = 2); col3 = Hausdorff
% outlier detect; col4 = RJE outlier detect (1 = valid, 0 = outlier)
%
% - will calculated PCI for L and R feet each as base, as well as for Merged series
% - all other time-domain cals will be performed by td_calcs
%
% The equations here assume that all outliers have already been marked (by
% td_outlier_detect). Much easier to do this first, because it requires both the
% time stamp channel and the voltage channel
%


if nargin < 2
    plot_phi = 0;
end

%% get the data set up

tM = Xmat(:,1); % all timestamps
chanM = Xmat(:,2); % channel

% get the final set of indices based on outlier detection, altrnating steps, longest run, etc.
valM = Xmat(:,end);


% do we have the correct number of channels?
num_chan = numel(unique(chanM));

if num_chan == 1
    % problem
    res.pciL = NaN;
    res.num_phiL = NaN;
    res.pciR = NaN;
    res.num_phiR = NaN;
    res.pciM = NaN;
    res.num_phiM = NaN;
    
else
    % OK
    
    tL = tM(chanM == 1);
    valL = valM(chanM == 1);

    tR = tM(chanM == 2);
    valR = valM(chanM == 2);



    %% Plotnik phi

    % rje decision: do phi on the ORIGINAL series, using NaN as a placeholder
    % for missing phi valMlues. THEN, delete the bad phi valMlues

    % to make life easier, let's re-assign L and R to "l" (longer) or "s"
    % (shorter), per the Plotnik paper

    % ***************
    % Left base
    tl = tL;
    vl = valL;
    ts = tR;
    vs = valR;

    tsmod = ts(ts > min(tl)); % only those ts events after the first tl
    minI = min(numel(tsmod),numel(tl));
    phi = nan(minI-1,1); % so we have NaN as placeholders, not 0

    l = 1; % need a separate counter for each foot
    s = 1;

    for i = 1:minI-1;

        if tl(l+1) > tsmod(s) && tsmod(s) > tl(l)
            % so far so good, but we also need to check if all steps are valMlid
            % (from outlier chanMlculation)

            check_out = [vl(i) vs(i) vl(i+1)];
            if sum(check_out) == 3
                % good
                phi(l+1) = 360 * ( tsmod(s) - tl(l) )/ ( tl(l+1) - tl(l) ); % per Plotnik 2007; the "+1" is so the index works correctly
            else
                phi(l+1) = NaN; % placeholder
            end

            l = l + 1;
            s = s + 1; 
        elseif tl(l) > tsmod(s)
            % this means we had two s's between a pair of Ls
            % need to keep the base the same and check the next opposite foot
            l = l + 0; % don't advalMnce to the next L pair, and don't put in NaN
            s = s + 1;
        elseif tsmod(s) > tl(l+1)
            % there is no s step here
            phi(l+1) = NaN; 
            l = l + 1;
            s = s + 0; % don't move
        end
    end

    phiL = phi; % store this

    % ***************
    % Right base

    % re-assign l and s
    tl = tR;
    vl = valR;
    ts = tL;
    vs = valL;


    tsmod = ts(ts > min(tl)); % only those ts events after the first tl
    minI = min(numel(tsmod),numel(tl));
    phi = nan(minI-1,1); % so we have NaN as placeholders, not 0

    l = 1; % need a separate counter for each foot
    s = 1;

    for i = 1:minI-1;

        if tl(l+1) > tsmod(s) && tsmod(s) > tl(l)
            % so far so good, but we also need to check if all steps are valMlid
            % (from outlier chanMlculation)

            check_out = [vl(i) vs(i) vl(i+1)];
            if sum(check_out) == 3
                % good
                phi(l+1) = 360 * ( tsmod(s) - tl(l) )/ ( tl(l+1) - tl(l) ); % per Plotnik 2007; the "+1" is so the index works correctly
            else
                phi(l+1) = NaN; % placeholder
            end

            l = l + 1;
            s = s + 1; 
        elseif tl(l) > tsmod(s)
            % this means we had two s's between a pair of Ls
            % need to keep the base the same and check the next opposite foot
            l = l + 0; % don't advance to the next L pair, and don't put in NaN
            s = s + 1;
        elseif tsmod(s) > tl(l+1)
            % there is no s step here
            phi(l+1) = NaN; 
            l = l + 1;
            s = s + 0; % don't move
        end
    end

    phiR = phi;

    % -----------------------------------------------------
    %% RJE phi method: a single series, ignoring gaps

    % merge tL and tR
    tM = [tL; tR];

    % ID L and R
    %iM = [repmat(1,size(tL)); repmat(2,size(tR))];

    % merge valL and valR
    valM = [valL; valR];

    % combine
    tv = [tM valM]; % 2 cols

    % re-sort
    tv = sortrows(tv,1);

    tM = tv(:,1);
    valM = tv(:,2);

    IntM = numel(tM); 
    phiM = nan(IntM-2,1);

    for m = 2:IntM-1;
       % some checks 
       valMl_check = [valM(m-1) valM(m-0) valM(m+1)];
       if sum(valMl_check) == 3; % if all index valMlues are good
           % ok
           phiM(m) = 360 * ( (tM(m) - tM(m-1))/ (tM(m+1) - tM(m-1))); 
       else
           phiM(m) = NaN;
           % and we advalMnce the counter one step automatichanMlly
       end
    end

    %% plot - with NaNs in place

    if plot_phi == 1;
        figure(100)

        subplot(3,1,1)
        plot(phiL)
        ylabel('phi (Left base)')
        xlabel('Ordinal Left events')

        subplot(3,1,2)
        plot(phiR)
        ylabel('phi (Right base)')
        xlabel('Ordinal Right events')

        subplot(3,1,3)
        plot(phiM)
        ylabel('phi (Merged)')
        xlabel('Ordinal Merged events')
    end


    %% chanMlculations: get rid of NaNs first
    % we can do this because each phi value is the result of the 3-number
    % calculation, so a valid phi means we have already ignored the bad points

    phiL(isnan(phiL)) = [];
    phiR(isnan(phiR)) = [];
    phiM(isnan(phiM)) = [];


    % next: other terms for PCI
    phiL_abs = mean(abs(phiL - 180));
    phiL_abs_prc = 100 * phiL_abs / 180;
    phiL_cv = 100 * std(phiL) / mean(phiL);

    phiR_abs = mean(abs(phiR - 180));
    phiR_abs_prc = 100 * phiR_abs / 180;
    phiR_cv = 100 * std(phiR) / mean(phiR);

    % RJE version of PCI
    %
    % 2013-08-11: let's do PCI as similar to Plotnik, but just using the whole
    % series. Note: we *cannot* do successive difference calculations, because
    % we have just removed all gaps in the phi series! Besides, phi is already
    % a successive difference measure. Save other successive difference
    % chanMlculations for PADM and PASD

    phiM_abs = mean(abs(phiM - 180));
    phiM_abs_prc = 100 * phiM_abs / 180;
    phiM_cv = 100 * std(phiM) / mean(phiM);

    % PCI
    pciL = phiL_abs_prc + phiL_cv; % this is on the scrubbed series
    pciR = phiR_abs_prc + phiR_cv;
    pciM = phiM_abs_prc + phiM_cv;


    %% outputs (valid)
    res.pciL = pciL;
    res.num_phiL = numel(phiL);
    res.pciR = pciR;
    res.num_phiR = numel(phiR);
    res.pciM = pciM;
    res.num_phiM = numel(phiM);
end
