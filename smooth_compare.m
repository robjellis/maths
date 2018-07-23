function output = smooth_compare(DTYPE,DPARAM,M,CV,N,SMTYPE,SMPARAM,ITER)

%
% smooth_compare(DTYPE,DPARAM,M,CV,N,SMTYPE,SMPARAM,ITER)
%
% generate noise (Gaussian, Poisson, or Gamma) and then smooth it (smooth_series.m) for as many
% values as are in SMPARAM. 


padl = nan(numel(N),numel(M),numel(CV),numel(SMPARAM),numel(ITER));
pasd = nan(numel(N),numel(M),numel(CV),numel(SMPARAM),numel(ITER));

if ITER > 10
   PLOT_IT = 0;
else
   PLOT_IT = 1;
end

j = 0; % for dfa counters

dfa_alpha = nan(ITER*numel(M)*numel(N)*numel(CV),numel(SMPARAM));

autocorr_vals = nan(ITER,20); % 

progressbar(0,0,0,0,0)

for n = 1:numel(N)
    progressbar((n-.0001)/numel(N),[],[],[],[])
    for m = 1:numel(M)
        progressbar([],(m-.0001)/numel(M),[],[],[])
        for c = 1:numel(CV)
            progressbar([],[],(c-.0001)/numel(CV),[],[])
            for s = 1:numel(SMPARAM)
                progressbar([],[],[],(s-.0001)/numel(SMPARAM),[])
                
                for i = 1:ITER
                    progressbar([],[],[],[],(i-.0001)/ITER)

                    res = smooth_series(DTYPE, DPARAM, M(m), CV(c),  SMTYPE, SMPARAM(s), N(n), PLOT_IT); 

                    Ygt = res.x_s;
                    Xgt = cumsum(Ygt);

                    % get the autocorrelation lags
                    vals = autocorr(Ygt);
                    vals = vals'; % row
                    vals = vals (2:end); % ignore lag0
                    
                    autocorr_vals(i,:) = vals;
                    
                    % stats for the *smoothed* series
                    do_dfa = 1;
                    plot_it = 0;
                    td_s = td_calcs_cum(Xgt,90,do_dfa,plot_it); 

                    % need to save a few things
                    dfa_alpha(i+j,s) = td_s.dfa_alpha;

                    padl(n,m,c,s,i) = td_s.padl;
                    pasd(n,m,c,s,i) = td_s.pasd;
                    
                end
                
                % for simplicity, just make a new figure every time we
                % finish ITER
                
                figure(350)
                subplot(2,3,s) % assumes 6 levels of SPMPARAM
                
                boxplot_rje(autocorr_vals,[2.5 97.5],0,0);
                
            end
            j = j + i; % constant incrementer
        end
    end
end


figure(250)

for s = 1:numel(SMPARAM)
    plot(padl(:),pasd(:),'.') % will turn them into a single row in the same way, so all values will correspond
    hold on
end
hold off

corr(padl(:),pasd(:),'type','spearman')

%% DFA boxplot
plot_outliers = 0;
output = boxplot_rje(dfa_alpha,[9 95],300,plot_outliers); % figure 300