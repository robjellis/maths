function output = mc_2samp(CO,EX,iter,nsamp,seed,log_trans,plot_figs)

% output = mc_2samp(CO,EX,iter,nsamp,seed,plot_figs)
%
% 2-sample tests (KS, AD, Mann Whitney U, U effect size) using N/2 random subsamples
%
% Note: it is assumed that CO and EX have the same exact dimensions; i.e.,
% column1 of CO will be compared with column1 of EX. This means that the
% code treats corresponding columns as *distinct measures* (e.g., age, IQ, reaction time, etc.)
%
% log_trans should have the same number of elements as there are columns in
% CO and EX
%
% http://www.mathworks.com/help/stats/ksdensity.html
%
% special case if 'iter' = 0; just take the original data without iterations
%
% updated = 2015.06.03 - option to log transform values prior to t-test 
% (1 = transform, 0 = don't transform)
%

tic;

% use a random seed
if nargin < 5 
    seed = randint(1,1,[1 1000000]);
end

% also ...
if seed == 0
    seed = randint(1,1,[1 1000000]);
end

if nargin < 6
    log_trans = zeros(size(EX,2));
end

if nargin < 7
    plot_figs = 0;
end

if iter == 0
    do_sub = 0;
    iter = 1;
elseif iter > 0
    do_sub = 1;
end


%% test options
do_rho = 0;
do_KS = 0;
do_AD = 0;

do_T = 1;
do_Z = 0; % don't worry about conversion to Z values
do_extra = 0; % extra stats from the U test

do_fwer = 0;


% get rid of NaN cases
co_nan = isnan(CO);
ex_nan = isnan(CO);
all_nan = sum(co_nan,2) + sum(ex_nan,2);
CO = CO(all_nan==0,:);
EX = EX(all_nan==0,:);

% new: change Inf to the highest observed value
   ex_vec = EX(:);
       maxd = max(ex_vec(ex_vec<Inf));
      isinf = EX == Inf;
EX(isinf) = maxd;

% size of D?
nsub = size(EX,1);
ncol = size(EX,2);

% make sure we have enough log_trans values
if max(log_trans) == 1
   if numel(log_trans) == ncol
       % ok
   else
       error(' Number of values in log_trans does not match number of measures.')
   end
end

% make sure same size as CO
if sum(size(CO) == size(EX)) == 2
    % ok
else
    error('Size of CO and EX do not match')
end



%% summary plot
%    if plot_figs == 1
%        figure(201)
%        for m = 1:ncol
%            d = data(:,m);
%            plot_this = 0;
%            fout = ecdf_manual(d,plot_this);
% 
%            if m == 1
%                 plot(fout.x,fout.f,'r','LineWidth',2)
%            else
%                 plot(fout.x,fout.f,'b')
%            end
%            hold on
%               
%        end
%        hold off
%        
%        xlabel('Data values')
%        ylabel('Cumulative Probability')
%    end


%% subsample analysis
% have to do this for each PAIR independently

KSsub = nan(iter,ncol); 
ADsub = nan(iter,ncol);
Tsub  = nan(iter,ncol);
Tsig  = nan(iter,ncol);

Usig  = nan(iter,ncol); % 1 if significant, 0 if not significant
Zsub  = nan(iter,ncol); % U to Z value
%Psub  = nan(iter,ncol); % p value for the approxmate (Z) test
PSsub = nan(iter,ncol); % the simple statistic: U / (n1*n2)
Asub  = nan(iter,ncol); % this is the "A" measure defined by Vargha and Delaney 2000 as a modification of the CLES
rho = nan(1,ncol);

% take the smaller value just in case we want to explore smaller subsamples    
n = min(nsamp,floor(nsub/2));

%progressbar(0,0) % *don't* initiliaze it

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% important: where do we want to take a new sample of subjects? Inside or
% outside the SF loop (columns)?

% set the random stream value; only have to do this once
% http://www.mathworks.com/help/matlab/ref/randstream.randperm.html
s = RandStream('mt19937ar','Seed',seed);

for i = 1:iter
    progressbar(i/iter) 
       
    if do_sub == 0
        i1 = 1:n; i1 = i1(:);
        i2 = n+1:n+n; i2 = i2(:);
    elseif do_sub == 1
        % permute the full set
        rp = randperm(s,nsub);
        
        % then take the first and second sets, which will be INDEPENDENT,
        % and will be used across the FULL set of columns
        i1 = rp(1:n); % for CO
        i2 = rp(n+1:n+n); % for EX
    end
    
    for m = 1:ncol
    %progressbar([],m/ncol) % takes too much time
    
        d1 = CO(:,m); 
        d2 = EX(:,m);
        
        % do we log transform this?
        log_this = log_trans(m);
        
        if log_this == 1
            % check for negative values
            if sum(CO < 0) > 0
                error(' Warning: negative values are present. Ln-transform not permitted.')
            end
        end
        
        
        if do_rho == 1
            % Spearman correlation for the full set of values
            if i == 1
                rho(m) = corr(d1,d2,'type','Spearman');      
            end
        end
        
        % get the subsample
        r1 = d1(i1);
        r2 = d2(i2);
        
    
        % KS test %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        

        if do_KS == 1
            [xxx xxx KSd] = kstest2(r1,r2);             
            KSsub(i,m) = KSd;
        end
        
        % AD test %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if do_AD == 1
            
            inds = [ones(numel(r1),1); ones(numel(r2),1)+1];
            r3 = [r1; r2];


            X = [r3 inds];

            % RJE modified AnDarksamtest to return "-1" if the algorithm fails
            % for some reason; this is useful because we can still tally these
            % values when we create the ECDF for this statistic!
        
            ADout = AnDarksamtest(X,.05);
            if ADout.ties_stat == -1
                ADsub(i,m) = NaN; % just exclude these cases
               % ADsub(i,m) = ADout.ties_stat;
            else
                % OK
                ADsub(i,m) = ADout.ties_stat;
            end
        end
        
        if do_T == 1
            % independent T test: ttest(B,A) will be positive if mean B > mean A
            if log_this == 0
                % don't do it
                [Th xxx xxx stats] = ttest2(r2,r1,.05,'both','unequal'); 
            elseif log_this == 1
                [Th xxx xxx stats] = ttest2(log(r2+1),log(r1+1),.05,'both','unequal'); 
            end
            %Tsub(i,m) = stats.tstat; % keep the sign so that we can tell about direction of inflation of Ex relative to Co!
            Tsub(i,m) = abs(stats.tstat); % abs to make things simpler
            Tsig(i,m) = Th; % 0 or 1
        end
        
        % Mann-Whitney U test: ranksum(B,A) will be positive Z if median B > median A
        % we actually use signrank_rje.m since it gives us several statistics we
        % need; i.e., we want to explicity quantify r2 (i.e., our Experimental value) with r1 (our control value)
        
        do_exact = 0; % don't get the actual exact P-value, just the exact binary significance (0 or 1)
        mw_out = ranksum_rje(r1,r2,do_exact);
        
        Usig(i,m)  = mw_out.H; % this is from exact test (hard coded by RJE)
        
        if do_extra == 1
            Zsub(i,m)  = abs(mw_out.Zex); % *not* the Matlab value of Z, but the standard U to Z conversion 
            %Psub(i,m)  = mw_out.Pexact; % two-tailed exact P
            PS = mw_out.PSex;

            % actually, let's do the "absolute" value of this
            PSsub(i,m) = 0.5 + abs(PS - 0.5);

            % the Vargha and Delaney 2000 effect size based upon Rex
            % RJE - let's come back to this later, may need to adjust!
            Rex = mw_out.Rex;
            A = (1/n) * (Rex/n - (n+1)/2);

            % but really let's make this "absolute value" to agree with the other measures!      
            Asub(i,m) = 0.5 + abs(A - 0.5);
        end

    end % column loop
    
end % iter loop

progressbar(1) % close it

% % we only care about the value of PS sub when the test is significant!
% 
% % get rid of the first column for simplicity
% PSsub(:,1) = NaN; % stay commented out!
% PSsub(Usig == 0) = NaN; 

% get the 95th percentile of PSsub
PSprc = prctile_nist(PSsub,95);

% one final thing: let's turn the Student T values into Z values so we
% don't have to worry about d.f.
df = 2 * n - 2;
pvals = tcdf(Tsub,df);
Tsub = norminv(pvals);

% we might have Inf here!
Tsub_vec = Tsub(:);

tm = max(Tsub_vec(Tsub_vec < Inf));

Tsub(Tsub == Inf) = tm;

if plot_figs == 1
    figure(202)
    clf

    figure(203)
    clf
end

if do_sub == 1
    %% get the observed p-values
    % let's make everything *positive* and one-tailed
    %prc1  = 95; % one tailed
    %prc2a = 2.5;
    %prc2b = 97.5; % two tailed

    % theoretical critical values
    KStheo = 1.36; % for > subjects per sample
    ADtheo = 2.492; % for >= 9 subjects per sample (per Pettitt 1976 Table 1, and Wiki)
    Ttheo = 1.96; % two-tailed alpha = .05; assuming we did T-to-Z
    Ztheo = 1.96; % two-tailed alpha = .05;
    Atheo = 0.64; % "medium effect" per Vargha Delaney Table 1
    
    % first we need the observed critical value for the first column
    if do_KS == 1
        KScrit = prctile_nist(KSsub(:,1),95);
    end
    
    if do_AD == 1
        ADcrit = prctile_nist(ADsub(:,1),95);
    end
    
    Tcrit_obs = prctile_nist(Tsub(:,1), 95);
    %Tlcrit = prctile_nist(Tsub(:,1), 2.5);
    %Thcrit = prctile_nist(Tsub(:,1),97.5);

    Zcrit = prctile_nist(Zsub(:,1), 95);
    %Zlcrit = prctile_nist(Zsub(:,1), 2.5);
    %Zhcrit = prctile_nist(Zsub(:,1),97.5);
    
    Acrit = prctile_nist(Asub(:,1), 95);

    KSsig = nan(1,ncol);
    ADsig = nan(1,ncol);
    %Tsig = nan(1,ncol);
     Zsig = nan(1,ncol);
     Asig = nan(1,ncol);
    

    % now look at the rest of the columns
    for m = 2:ncol
        
        ks_cur = KSsub(:,m);
        ad_cur = ADsub(:,m);
         t_cur =  Tsub(:,m);
         z_cur =  Zsub(:,m);
         a_cur =  Asub(:,m);
        ps_cur = PSsub(:,m);

        % get the upper tail (for plotting purposes)
        %ks_prc = prctile_nist(ks_cur,95);
        %ad_prc = prctile_nist(ad_cur,95);
        
        % need lower and upper tails for plotting
        %t_prc   = prctile_nist(t_cur, 95);
        %tl_prc  = prctile_nist(t_cur, 2.5);
        %th_prc  = prctile_nist(t_cur,97.5);
        
        %z_prc   = prctile_nist(z_cur, 95);
        %zl_prc  = prctile_nist(z_cur, 2.5);
        %zh_prc  = prctile_nist(z_cur,97.5);
        
        % now tally how many results are significant relative to the observed critical value
        KSsig(m) = 100 * sum(ks_cur > KStheo) / numel(ks_cur);
        ADsig(m) = 100 * sum(ad_cur > ADtheo) / numel(ad_cur);
        
        % new - just directly take the significance from the ttest!
        %Tsig(m) = 100 * ( sum(t_cur > Thcrit) + sum(t_cur < Tlcrit) )  / numel(t_cur);
        %Tsig(m) = 100 * ( sum(t_cur > Ttheo) ) / numel(t_cur);
         
        %Zsig(m) = 100 * ( sum(z_cur > Zhcrit) + sum(t_cur < Zlcrit) )  / numel(z_cur);
         Zsig(m) = 100 * ( sum(z_cur > Ztheo) ) / numel(z_cur);
         
         Asig(m) = 100 * ( sum(a_cur > Atheo) ) / numel(a_cur);

        %% figure 202 - don't need these
%         figure(202)
%         subplot(2,1,1)
%         if m == 2 % plot this last
%             eo = ecdf_manual(KSsub(:,1),0);
%             plot(eo.x,eo.f,'r','LineWidth',3)
%         end
%         hold on
%             ecdf_manual(ks_cur,1);        
%             %plot([ks_prc ks_prc],[0 1],'g','LineWidth',2)
%             %plot([KScrit KScrit],[0 1],'r')  
%         hold off
%         ylabel('KS')
% 
%         subplot(2,1,2)
%         hold on
%             if m == 2 % plot this last
%                 eo = ecdf_manual(Tsub(:,1),0);
%                 plot(eo.x,eo.f,'r','LineWidth',3)
%             end
%             ecdf_manual(t_cur,1);
%             %plot([t_prc t_prc],[0 1],'g','LineWidth',2)
%             %plot([tl_prc tl_prc],[0 1],'g','LineWidth',2)
%             %plot([th_prc th_prc],[0 1],'g','LineWidth',2)
%             %plot([Tcrit_obs Tcrit_obs],[0 1],'r')
%             %plot([Tlcrit Tlcrit],[0 1],'r') 
%             %plot([Thcrit Thcrit],[0 1],'r') 
%         hold off    
%         %axis([0 tmax 0 1])
%         ylabel('T-to-Z')
        

        %% figure 203
        if plot_figs == 1
        figure(203)
        
        subplot(1,4,1)
        hold on
            if m == 2 % plot this last
                %eo = ecdf_manual(ADsub(:,1),0);
                %plot(eo.x,eo.f,'r','LineWidth',3) % plot the "null" curve
            end        
            ecdf_manual(ad_cur,1);
            %plot([ad_prc ad_prc],[0 1],'g','LineWidth',2)
            %plot([ADcrit ADcrit],[0 1],'r') 
        hold off
        xlabel('AD statistic')
        ylabel('Cum. Prob.')
        % don't set an axis limit in case we have "-1" value from the test

        
        subplot(1,4,2)
        hold on
            if m == 2 % plot this last
                eo = ecdf_manual(Zsub(:,1),0);
                plot(eo.x,eo.f,'r','LineWidth',3) % plot the "null" curve
            end
            ecdf_manual(z_cur,1);
            %plot([z_prc z_prc],[0 1],'g','LineWidth',2)
            %plot([zl_prc zl_prc],[0 1],'g','LineWidth',2)
            %plot([zh_prc zh_prc],[0 1],'g','LineWidth',2)
            %plot([Zcrit Zcrit],[0 1],'r') 
            %plot([Zlcrit Zlcrit],[0 1],'r') 
            %plot([Zhcrit Zhcrit],[0 1],'r') 
        hold off
        axis([0 max(Zsub(:)) 0 1])
        
        xlabel('U-to-Z statistic')
        ylabel('Cum. Prob.')
        
        
        subplot(1,4,3)
        hold on
            if m == 2 % plot this last
                eo = ecdf_manual(Asub(:,1),0);
                plot(eo.x,eo.f,'r','LineWidth',3) % plot the "null" curve
            end
            ecdf_manual(a_cur,1,'b'); 
            %plot([z_prc z_prc],[0 1],'g','LineWidth',2)
            %plot([Zcrit Zcrit],[0 1],'r') 
            %plot([0.5 0.5], [0 1], 'c') % just to highlight the neutral point
        hold off
        axis([0.5 1.0 0 1]) % with the "absolute value" of A, we never go below 0.5
        xlabel('A (effect size) statistic')       
        ylabel('Cum. Prob.')
        
        
        subplot(1,4,4)
        hold on
            if m == 2 % plot this last
                eo = ecdf_manual(PSsub(:,1),0);
                plot(eo.x,eo.f,'r','LineWidth',3) % plot the "null" curve
            end
            ecdf_manual(ps_cur,1,'b'); 
            %plot([z_prc z_prc],[0 1],'g','LineWidth',2)
            %plot([Zcrit Zcrit],[0 1],'r') 
            %plot([0.5 0.5], [0 1], 'c') % just to highlight the neutral point
        hold off
        axis([0.5 1.0 0 1]) % with the "absolute value" of A, we never go below 0.5
        xlabel('PS statistic')       
        ylabel('Cum. Prob.')        
        %pause(.1)
        
        %% figure 204 - correlations
        
        figure(204)
        % PS vs Z
        subplot(1,2,1)
        plot(PSsub(:),Zsub(:),'.')
        
        % PS vs A
        subplot(1,2,2)
        plot(PSsub(:),Asub(:),'.')
        
        end
        
    end

    %% 2014.05.26 - FWER test
    % 1. rank order all measures by FPR, and SORT the Usig matrix 
    % 2. for n = 1:10, just take max(rows) to get the FWER
    % 3. do the same for the HIGEST FPR values
    
    % This is SPECIFICALLY designed for the HRV project
    
    %Hsig = sum(Usig) / iter;
    
    %[Hsort Hind] = sort(Hsig);
    
    %Usig_sort = Usig(:,Hind); % now this is sorted
    
    if do_fwer == 1
        if ncol == 33
            cols = [2 4 23 26 28 30 29 31];
            FWER_t = nan(4,1);
            FWER_u = nan(4,1);

            for r = 2:2:numel(cols);
                Tdata = Tsig(:,cols(1:r));
                Udata = Usig(:,cols(1:r));

                FWER_t(r/2) = sum(max(Tdata,[],2)) / iter; 
                FWER_u(r/2) = sum(max(Udata,[],2)) / iter; 
            end

        else
            FWER_t = NaN;
            FWER_u = NaN;
        end
    
    end

%% agreement of significance tests
% rje: in theory, we could do the Pearson Phi coefficient, but a problem
% arises when we have all 0s or 1s (which can happen with small number of
% iterations

% a simpler way is just to get the proportion of iterations in which the
% significance of x and y matches!


%% 2014.05.30 - the Phi coefficient
% just a Pearson correlation on binary values across the set of HRV measures

if ncol == 33
    use_col = [1:4 7:ncol]; % ignore PNN50 and PNN40
else
     use_col = 1:ncol;
end
  
ncol = numel(use_col);
Tuse = Tsig(:,use_col);
Uuse = Usig(:,use_col);

combs = combnk(1:ncol,2);
ncombs = size(combs,1);

if do_extra == 1
    Tmatch = nan(ncombs,1);
    Umatch = nan(ncombs,1);

    for c = 1:ncombs

        Tx = Tuse(:,combs(c,1));
        Ty = Tuse(:,combs(c,2));

        Ux = Uuse(:,combs(c,1));
        Uy = Uuse(:,combs(c,2));

        Tmatch(c) = sum(Tx == Ty) / iter;
        Umatch(c) = sum(Ux == Uy) / iter;

    end  

    % 2014.06.02 - putting corr and cov back in the game
    % this SHOULD work so long as we do at least 50 iterations so that we
    % observe at least some false positives

    corr_t = corr(Tsig);
    corr_u = corr(Usig);

    cov_t = cov(Tsig);
    cov_u = cov(Usig);
end


    % there are a few special cases in which doing corr(x,y) won't work
    % (nor does the manual formula for Pearson phi)
    
%        % for T
%        if sum(x==0) == iter && sum(y==0) == iter % if both are all 0s
%             phi = 1;
%        elseif sum(x==1) == iter && sum(y==1) == iter % if both are all 1s
%             phi = 1;
%        elseif sum(x==1) == iter && sum(y==0) == iter
%             phi = -1;
%        elseif sum(x==0) == iter && sum(y==1) == iter 
%             phi = -1;
%        else
%            % in theory, if either x or y is constant, then we will still
%            % get an NaN, but this is unlikely; all 0s in either vector is the most likely
%            % scenario
%            phi = corr(x,y);
%        end
       
%        if j == 1
%            phi_t(c) = phi;
%        elseif j == 2
%            phi_u(c) = phi;
%        end
       
       

% % get the lower half, and EXCLUDE the main diagonal (r = 1.0000)
% phi_t = tril(phi_t,-1);
% phi_u = tril(phi_u,-1);
% 
% % vector
% phi_t = phi_t(:);
% phi_u = phi_u(:);
% 
% % get rid of zeros
% 
% phi_t(phi_t == 0) = [];
% phi_u(phi_u == 0) = [];


    
    %% output

    output.n      = n;
    output.iter   = iter;
    
    if do_KS == 1
        output.KStheo = KStheo;
        output.KScrit = KScrit; % just to compare with the theoretical value
        output.KSsig = KSsig;
    end
    
    if do_AD == 1
        output.ADtheo = ADtheo;
        output.ADcrit = ADcrit; % just to compare with the theoretical value
        output.ADsig = ADsig;
    end
    
    %output.Ttheo = Ttheo;
    %output.Tcrit_obs = Tcrit_obs; % just to compare with the theoretical value
    output.Tsig_all = Tsig; % the original data
    output.Tsig_FPR  = sum(Tsig) / iter; % ranges from 0 to 1

    if do_Z == 1
        output.Ztheo = Ztheo;
        output.Zcrit = Zcrit; % just to compare with the theoretical value
        output.Zsig = Zsig;
    end

       
    if do_extra == 1
        output.Atheo = Atheo;
        output.Asig = Asig;

        output.PSsub = PSsub; % the whole matrix
        output.PSprc = PSprc; % 95th percentile
    end
    
    % do this as a PROPORTION
    output.Usig_all = Usig;
    output.Usig_FPR = sum(Usig) / iter; % MWW test significance
    
    if do_fwer == 1
        output.FWER_t = FWER_t;
        output.FWER_u = FWER_u;
    end
    
    if do_rho == 1
        output.rho = rho;
    end
    
    if do_extra == 1
        % proportion of iterations that have same sig (0 or 1)
        output.Tmatch = Tmatch;
        output.Umatch = Umatch;

        % corr and cov
        output.corr_t = corr_t;
        output.corr_u = corr_u;
        output.cov_t = cov_t;
        output.cov_u = cov_u;
    end
    

else
    
    output.AD = ADsub;
    output.T = Tsub;
    output.U = Usub;
    output.Uz = Zsub;
    output.A = Asub;
    
    
end % do_sub 

ttoc = toc;
output.time = ttoc;

