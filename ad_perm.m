function output = ad_perm(data,iter,plot_figs,label)

% 2-sample KS and AD test with permutation testing for critical values

% http://www.mathworks.com/help/stats/ksdensity.html

if nargin < 2
   iter = 1000;
   plot_figs = 1;
   label = '';
end

% get rid of NaN cases
data_nan = isnan(data);
data_nan = sum(data_nan,2);

data = data(data_nan==0,:);
data_vec = data(:);

% size of D?
nsub = size(data,1);
ncol = size(data,2);

% get the full range of the data
xmin = min(data_vec);
xmax = max(data_vec);
range = xmax - xmin;


%% calculate the vertical distance D for each pair


KSstat = nan(ncol,1);
KSpval = nan(ncol,1);

ADstat = nan(ncol,1);
ADpval = nan(ncol,1);

% get the unique pairs
pairs = nchoosek(1:ncol,2);

npairs = size(pairs,1);

for p = 2:ncol
   %d1 = data(:,pairs(p,1));
   %d2 = data(:,pairs(p,2));
   d1 = data(:,1); % always the reference data
   d2 = data(:,p);
   
   % alternative: use RJE function ecdf_manual

   [xxx KSp KSd] = kstest2(d1,d2);
   
   %KSstat(pairs(p,2),pairs(p,1)) = KSd;
   %KSpval(pairs(p,2),pairs(p,1)) = KSp;
   
   KSstat(p) = KSd;
   
   inds = [ones(numel(d1),1); ones(numel(d2),1)+1];
   d3 = [d1; d2];
   X = [d3 inds];
   
   % AD test, adjusted for ties
   ADout = AnDarksamtest(X,.05);
   %ADstat(pairs(p,2),pairs(p,1)) = ADout.ties_stat;
   %ADpval(pairs(p,2),pairs(p,1)) = ADout.ties_pval;
   
   ADstat(p) = ADout.ties_stat;
   
   if plot_figs == 1
       figure(201)
       plot_this = 0;
       f1out = ecdf_manual(d1,plot_this);
       f2out = ecdf_manual(d2,plot_this);
       
       plot(f1out.x,f1out.f,'r','LineWidth',2)
       hold on
       plot(f2out.x,f2out.f,'b')
       
       xlabel('Data values')
       ylabel('Cumulative Probability')
   end

end % column loop
hold off

%% now we need to find the critical value using bootstrap analysis
% have to do this for each PAIR independently (permutation testing)

KSperm = nan(iter,ncol); 
ADperm = nan(iter,ncol);

progressbar(0,0) % initialize

    % this will be for the permuted data
    p1 = nan(nsub,1);
    p2 = nan(nsub,1);
    
for i = 1:iter
    progressbar([],i/iter)
    
    c = randint(nsub,1,[1 2]); % unique for each iter
  
    for m = 2:ncol
    progressbar((m-1)/ncol,[])
     
        % pull a single value randomly from each *row*
        %c1 = randint(nsub,1,[1 ncol]);
        %c2 = randint(nsub,1,[1 ncol]);

        % get the real data        
        r = [data(:,1) data(:,m)];

        % this will be for the permuted data
        p1 = nan(nsub,1);
        p2 = nan(nsub,1);
    
        % get the two vectors
        for s = 1:nsub
            
            if c(s) == 1
                p1(s) = r(s,1);
                p2(s) = r(s,2);
            elseif c(s) == 2
                p1(s) = r(s,2);
                p2(s) = r(s,1);                
            end
            
        end

        % KS test
        [xxx xxx KSd] = kstest2(p1,p2);
              
        KSperm(i,m) = KSd;
        
        % AD test
        inds = [ones(numel(p1),1); ones(numel(p2),1)+1];
        p3 = [p1; p2];
        X = [p3 inds];
   

        ADout = AnDarksamtest(X,.05);
        ADperm(i,m) = ADout.ties_stat;

    end % column loop
    
end % iter loop
progressbar(1)

%% get the observed p-values

for m = 2:ncol
    
    ks_cur = KSperm(:,m);
    ad_cur = ADperm(:,m);
    
    % get the upper tail
    ks_prc = prctile_nist(ks_cur,95);
    ad_prc = prctile_nist(ad_cur,95);
    
    KSpval(m) = sum(ks_cur >= KSstat(m,1)) / numel(ks_cur);
    ADpval(m) = sum(ad_cur >= ADstat(m,1)) / numel(ad_cur);
    
%     % now get the PDF    
%     f_crit  = ksdensity(ovals_cur,xi_crit,'width',bw_crit_min);
% 
%     % now get the CDF
%     cdf_crit  = ksdensity(ovals_cur,xi_crit,'width',bw_crit_min,'function','cdf');
% 
%     % and get the one-tailed critical value (p = .05)
%     x_crit_prc(m) = prctile_nist(ovals_cur,5);
%     x_crit_kde(m) = max(xi_crit(cdf_crit <= .05));

    figure(203)
    subplot(2,1,1)
    ecdf_manual(ks_cur,1)
    hold on
        plot([ks_prc ks_prc],[0 1],'g','LineWidth',2)
        plot([KSstat(m,1) KSstat(m,1)],[0 1],'r')  
    hold off
    
%     hmax = hist(ovals_cur,20);
%     hmax = max(hmax);
%     
%     hold on
%     plot([x_crit_prc(m) x_crit_prc(m)],[0 10],'r')
%     set(gca,'xlim',[xmin xmax])
%     hold off
    
    subplot(2,1,2)
    ecdf_manual(ad_cur,1)
    hold on
        plot([ad_prc ad_prc],[0 1],'g','LineWidth',2)
        plot([ADstat(m,1) ADstat(m,1)],[0 1],'r') 
    hold off
%     plot(xi_crit,f_crit,'LineWidth',2)
%     hold on
%     plot([x_crit_kde(m) x_crit_kde(m)],[0 max(f_crit)],'r')
%     hold off
    pause(1)
end


%% output

output.KSstat = KSstat;
output.KSpval = KSpval;
output.ADstat = ADstat;
output.ADpval = ADpval;

% output.data_bw_kde = bw_kde;
% output.data_bw_matlab = bw_matlab;
% output.data_bw_final = bw_final;
% output.data_ovl = D_mat;
% 
% output.perm_bw = bw_crit_min;
% output.perm_ovl = ovals;
% 
% output.x_crit_kde = x_crit_kde;
% output.x_crit_prc = x_crit_prc;
% output.ovl_gt = D_mat(:,1); % versus ground truth (column 1)
% output.sig_kde = output.ovl_gt < x_crit_kde;
% output.sig_prc = output.ovl_gt < x_crit_prc;




