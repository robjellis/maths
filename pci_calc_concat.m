function [res] = pci_calc_concat(Xmat,plot_phi, phi_max_ev, min_L_R, toe_off)

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
% * Note: An EVEN number of total events is recommended, so that we
% maximize the number of PHI values that can be calculated (even means L
% and R phi will utilize the full data set)
%
% Xmat is a four-column matrix that is output from td_outlier_detect by
% RJE. col1 = *timestamps*; col2 = channel (L = 1, R = 2); col3 = Hausdorff
% outlier detect; col4 = RJE outlier detect (1 = valid, 0 = outlier)
% 
% * revised: Xmat can be a single column; we'll just take alternating
% events as Right (odd) and Left (even)
%
% - will calculated PCI for L and R feet each as base, as well as for Merged series
% - all other time-domain cals will be performed by td_calcs
%
% * phi_max_ev will take the first N events and perform PCI calc just on
% that data; if phi_max_ev = Inf, then all data will be used
%
% *** The equations here assume that all outliers have already been marked (by
% td_outlier_detect). Much easier to do this first, because it requires both the
% time stamp channel and the voltage channel ***
%
%
% rje | 2013.09.15
% update 2014.03.06: now able to be called by td_calcs_concat


if nargin < 2
    plot_phi = 0;
end

if nargin < 3
    phi_max_ev = Inf;
end

if nargin < 4 || isempty(toe_off)
    do_swing = 0;
else
    do_swing = 1;
end

%% make sure we have two channels
% revised 2014.03.06: if we just have one channel we take odd as Right and even as Left

Xmat_cols = size(Xmat,2);

if Xmat_cols == 4 && numel(unique(chan_check)) == 2
    % OK
elseif Xmat_cols == 1
   % create the matrix; we assume ALL DATA points are valid
   Xmat_new = nan(numel(Xmat),4);
   Xmat_new(:,1) = Xmat;
   Xmat_new(:,3) = ones(numel(Xmat,1));
   Xmat_new(:,4) = ones(numel(Xmat,1));
   
   % create the alternating series
   alt = zeros(numel(Xmat),1);
   alt(1:2:end) = 2; % right
   alt(2:2:end) = 1; % left
   
   Xmat_new(:,2) = alt;
   
   Xmat = Xmat_new;
end


%% get the locations of the run(s) (always from the last column)

% note: if we input a *single* vector of numbers, it is assumed that ALL
% data points are valid (i.e., we just extract the longest stable segment
% for the BEATS project). So, we need a little help here ...

% if min(size(Xmat)) == 1;
%    run_inds = ones(size(Xmat,1),1);
% else
%    run_inds = Xmat(:,end);
% end

run_inds = Xmat(:,end);

run_inds = [0; run_inds; 0]; % the zeros allow us to identify beginning ([0 1 1 1]) or ending (e.g., [1 1 1 0]) runs
diffrun = diff(run_inds); % indexes will be correct to mark timestamps
    
startIndex = find(diffrun == 1); % find ALL changes from 0 to 1

endIndex = find(diffrun == -1); % final ALL changes from 1 to 0, or a 1 at end of vector
        
endIndex = endIndex - 1; % the last 1 before a 0; the "-1" makes this work correctly since we are using diff()

nruns = numel(startIndex); % so we can loop through this

% set up placeholders
phiL = [];
phiR = [];
phiM = nan(size(Xmat,1),1);
swingL = [];
swingR = [];
swingM = []; % just so we capture this in correct temporal order

% get the toe off times
if do_swing == 1
    toeL = toe_off.L; toeL = toeL(:);
    toeR = toe_off.R; toeR = toeR(:);
    toeM = [toeL; toeR]; toeM = sort(toeM);
end

%% loop
for n = 1:nruns
    sL = []; % just to clear it
    sR = [];
    sM = [];
    
    % get the "slice" of Xmat we need; note: we have already done outlier
    % detection, altrnating steps, longest run, etc.
    tM    = Xmat(startIndex(n):endIndex(n),1); % all timestamps
    chanM = Xmat(startIndex(n):endIndex(n),2); % channel

    %valM = Xmat(:,end); % get the final set of indices

    
    tL = tM(chanM == 1);
    %valL = valM(chanM == 1);

    tR = tM(chanM == 2);
    %valR = valM(chanM == 2);    
    
    % determine if we do phi calculations
    if numel(tL) >= min_L_R && numel(tR) >= min_L_R; % need at least this many events in a consecutive string
        do_phi = 1;
    else
        do_phi = 0;
    end
    
    if do_swing == 1
        %% First, we calculate swing times
        for i = 1:numel(tL)-1
            ind = min(intersect(find(toeL > tL(i)),find(toeL < tL(i+1))));

            if isempty(ind)
                sL(i) = NaN; % placeholder
            else
                sL(i) = tL(i+1) - toeL(ind);
            end
        end
        sL = sL(:);

        for i = 1:numel(tR)-1
            ind = min(intersect(find(toeR > tR(i)),find(toeR < tR(i+1))));

            if isempty(ind)
                sR(i) = NaN; % placeholder
            else
                sR(i) = tR(i+1) - toeR(ind);
            end        
        end
        sR = sR(:);

        % merged, in case we want this
        % (since we've already made sure that we have L and R alternating, this method is fine)

        % first pair the timestamp with the swing time
        Lpair = [tL(2:end) sL];
        Rpair = [tR(2:end) sR];

        % combine them
        Mpair = [Lpair; Rpair];
        Mpair = sortrows(Mpair,1); % will sort by the timestamp column

        sM = Mpair(:,2);

        % concatenate across runs
        swingL = [swingL; sL];
        swingR = [swingR; sR];
        swingM = [swingM; sM];
    end % do_swing == 1

    %% Plotnik phi
    if do_phi == 1
    % rje decision: do phi on the ORIGINAL series, using NaN as a placeholder
    % for missing phi values. THEN, delete the bad phi values

    % to make life easier, let's re-assign L and R to "l" (longer) or "s"
    % (shorter), per the Plotnik paper

    % ***************
    % Left base
    tl = tL;
    ts = tR;

    % only those ts events after the first tl and before the last tl
    tsmod = ts(ts > min(tl)); 
    tsmod = tsmod(tsmod < max(tl));
    
    if numel(tsmod) == numel(tl)
        tsmod = tsmod(1:end-1);
        warning('fix me!')
    end
    minI = numel(tsmod);
    
    % reset to be safe
    phi = [];

    for i = 1:minI;

        phi(i,1) = 360 * ( tsmod(i) - tl(i) )/ ( tl(i+1) - tl(i) ); % per Plotnik 2007
 
    end
    
    phiL = [phiL; phi; nan(2,1)]; % store this; the NaNs are just so we can see the GAPS clearly

    % ***************
    % Right base

    % re-assign l and s
    tl = tR;
    ts = tL;

    % only those ts events after the first tl and before the last tl
    tsmod = ts(ts > min(tl)); 
    tsmod = tsmod(tsmod < max(tl));
    minI = numel(tsmod);

    % reset to be safe
    phi = [];

     for i = 1:minI;

         phi(i,1) = 360 * ( tsmod(i) - tl(i) )/ ( tl(i+1) - tl(i) ); % per Plotnik 2007

     end

    phiR = [phiR; phi; nan(2,1)]; % store this; the NaNs are just so we can see the GAPS clearly

    
    % -----------------------------------------------------
    % RJE phi method: a single series, ignoring gaps

    IntM = numel(tM); 

    for m = 2:IntM-1;

       phiM(startIndex(n)+m-1) = 360 * ( (tM(m) - tM(m-1))/ (tM(m+1) - tM(m-1))); 

    end
    
    else
        % just add a gap
        phiR = [phiR; nan(2,1)]; % the NaNs are just so we can see the GAPS clearly
        phiL = [phiL; nan(2,1)]; % the NaNs are just so we can see the GAPS clearly
        phiM = [phiM; nan(2,1)]; % the NaNs are just so we can see the GAPS clearly
    end

end % run loop

% 2014.03.12
% need to get rid of NaNs!
phiL = phiL(isnan(phiL)==0);
phiR = phiR(isnan(phiR)==0);
phiM = phiM(isnan(phiM)== 0);

% for simplicity sake, we want to have same N in phiL and phiR
nmin = min(numel(phiL),numel(phiR));

phiL = phiL(1:nmin);
phiR = phiR(1:nmin);
   
    %% plot - with NaNs in place

    if do_phi == 1 && plot_phi == 1;
        figure(100)

        subplot(2,1,1)
        plot(phiL,'Marker','.','Color','r') % we color this in red, since R is the *variable* and L is the *base*
        hold on
        plot(phiR,'Marker','.','Color','b')
        xlabel('Ordinal events')
        ylabel('phi (red = right base; blue = left base)')
        hold off
        
        subplot(2,1,2)
        plot(phiM,'Marker','.','Color',[0 0 0])
        ylabel('phi (Merged)')
        xlabel('Ordinal Merged events')
        
        min_ev = min(numel(phiL),numel(phiR));
        figure(101)
        plot(phiL(1:min_ev),phiR(1:min_ev),'.')
        
    end


    %% calculations: get rid of NaNs first
    % we can do this because each phi value is the result of the 3-number
    % calculation, so a valid phi means we have already ignored the bad points

    phiL(isnan(phiL)) = [];
    phiR(isnan(phiR)) = [];
    phiM(isnan(phiM)) = [];
    
    if do_swing == 1
        swingL(isnan(swingL)) = [];
        swingR(isnan(swingR)) = [];
        swingM(isnan(swingM)) = [];


        plot_swing = 0;
        if plot_swing == 1
            figure(410)
            subplot(3,1,1)
            plot(swingL)
            subplot(3,1,2)
            plot(swingR)
            subplot(3,1,3)
            plot(swingM)
        end
        
    else
        swingL = [];
        swingR = [];
        swingM = [];
    end
    
    %% now we cap the number of events
    
    % first, how many events do we actually have?  
    res.phiL_num = numel(phiL);
    res.phiR_num = numel(phiR);
    res.phiM_num = numel(phiM);
    
    num_events = phi_max_ev * 2;
    
    if phi_max_ev == Inf
        % keep all events; don't need to redefine
   
    else
        if phi_max_ev <= numel(phiL)
            phiL = phiL(1:phi_max_ev);
        %else
            %phiL = NaN;
        end
        
        if phi_max_ev <= numel(phiR)
            phiR = phiR(1:phi_max_ev);
        %else
            %phiR = NaN;
        end
        
        if num_events <= numel(phiM)
            phiM   = phiM(1:num_events);  
        %else
            %phiM = NaN;
        end
        
        if phi_max_ev <= numel(swingL)
            swingL   = swingL(1:phi_max_ev);  
        %else
            %swingL = NaN;
        end
        
        if phi_max_ev <= numel(swingR)
            swingR   = swingR(1:phi_max_ev);  
        %else
            %swingR = NaN;
        end
        
        if num_events <= numel(swingM)
            swingM   = swingM(1:num_events);  
        %else
            %swingM = NaN;
        end        
    end
    
    
    % next: other terms for PCI
    phiL_sd = std(phiL);
    phiL_mn = mean(phiL);
    phiL_abs = mean(abs(phiL - 180));
    phiL_abs_prc = 100 * phiL_abs / 180;
    phiL_cv = 100 * phiL_sd / phiL_mn;
    
    phiR_sd = std(phiR);
    phiR_mn = mean(phiR);
    phiR_abs = mean(abs(phiR - 180));
    phiR_abs_prc = 100 * phiR_abs / 180;
    phiR_cv = 100 * phiR_sd / phiR_mn;

    %% RJE version of PCI
    %
    % 2013-08-11: let's do PCI as similar to Plotnik, but just using the whole
    % series. Note: we *cannot* do successive difference calculations, because
    % we have just removed all gaps in the phi series! Besides, phi is already
    % a successive difference measure. Save other successive difference
    % calculations for PADM and PASD

    phiM_sd = std(phiM);
    phiM_mn = mean(phiM);
    phiM_abs = mean(abs(phiM - 180));
    phiM_abs_prc = 100 * phiM_abs / 180;
    phiM_cv = 100 * phiM_sd / phiM_mn;

    % PCI
    pciL = phiL_abs_prc + phiL_cv; % this is on the scrubbed series
    pciR = phiR_abs_prc + phiR_cv;
    pciM = phiM_abs_prc + phiM_cv;

    
    %% now simple calculations for swing time
    swingL_mn = mean(swingL);
    %swingL_sd = std(swingL);
    swingL_cv = 100 * std(swingL) / mean(swingL);
    swingR_mn = mean(swingR);
    %swingR_sd = std(swingR);
    swingR_cv = 100 * std(swingR) / mean(swingR);
    swingMer_mn = mean(swingM);
    %swingMer_sd = std(swingM);
    swingMer_cv = 100 * std(swingM) / mean(swingM);
    
    swingMax_mn = max(swingL_mn,swingR_mn);
    %swingMax_sd = max(swingL_sd,swingR_sd);
    swingMax_cv = max(swingL_cv,swingR_cv);
    
    %% outputs (valid; NaNs defined above)
    res.numL = numel(tL);
    res.numR = numel(tR);
    res.tL = tL;
    res.tR = tR;

    res.cvL = 100* std(diff(tL)) / mean(diff(tL));
    res.cvR = 100* std(diff(tR)) / mean(diff(tR));
    
    rmsL = rmssd(diff(tL));
    rmsR = rmssd(diff(tR));
    res.cvrmsL = rmsL.cvrms;
    res.cvrmsR = rmsR.cvrms;
    
    res.phiL = phiL; % the full series
    res.phiR = phiR;
    res.phiM = phiM;

    res.phiL_mn = phiL_mn;
    res.phiR_mn = phiR_mn;
    res.phiM_mn = phiM_mn;
    
    res.phiL_abs_prc = phiL_abs_prc;
    res.phiR_abs_prc = phiR_abs_prc;
    res.phiM_abs_prc = phiM_abs_prc;
    
    res.phiL_sd = phiL_sd;
    res.phiR_sd = phiR_sd;
    res.phiM_sd = phiM_sd;
    
    res.phiL_cv = phiL_cv;
    res.phiR_cv = phiR_cv;
    res.phiM_cv = phiM_cv;
    
    % RJE version which assumes mn = 180
    res.phiL_altcv = 100* phiL_sd / 180;
    res.phiR_altcv = 100* phiR_sd / 180;
    res.phiM_altcv = 100* phiM_sd / 180;
    
    res.pciL = pciL;
    res.num_phiL = numel(phiL);
    res.pciR = pciR;
    res.num_phiR = numel(phiR);
    res.pciMer = pciM;
    res.num_phiMer = numel(phiM);
    
    res.swingL_num = numel(swingL_mn);
    res.swingL_mn = swingL_mn;
    %res.swingL_sd = swingL_sd;
    res.swingL_cv = swingL_cv;
    
    res.swingR_num = numel(swingR_mn);
    res.swingR_mn = swingR_mn;
    %res.swingR_sd = swingR_sd;
    res.swingR_cv = swingR_cv;
    res.swingMer_mn = swingMer_mn;
    %res.swingMer_sd = swingMer_sd;
    res.swingMer_cv = swingMer_cv;
    res.swingMax_mn = swingMax_mn;
    %res.swingMax_sd = swingMax_sd;
    res.swingMax_cv = swingMax_cv;

