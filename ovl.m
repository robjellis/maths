function output = ovl(data,iter,plot_figs)

% overlap coefficients

% http://www.mathworks.com/help/stats/ksdensity.html

% get rid of NaN cases
data_nan = isnan(data);
data_nan = sum(data_nan,2);

data = data(data_nan==0,:);
data_vec = data(:);

% size of D?
nsub = size(data,1);
ncol = size(data,2);

% get min and max of dataset
% MIN=min(data)-Range/10 and MAX=max(data)+Range/10; per Botev code
xmin = min(data_vec);
xmax = max(data_vec);
range = xmax - xmin;

range_div = 4;

% we CAN'T do the procedure below, because if we do, then the inidividual
% distributions won't have an area of 1.0, which makes the intersection not
% valid!
% % restrict the range if needed for a few special cases
% 
% if xmin >= 0
%     % thre are no negative numbers (like RMSSD) or we might have a percent
%     % score, so we don't want to go negative
%     xmin = max(0, xmin - range/range_div); % take the larger of the two
% end
% 
% if xmax <= 100
%     % might be a percent score
%     xmax = min(100, xmax + range/range_div); % take the lesser value
% end

xmin = xmin - range/range_div;
xmax = xmax + range/range_div;

% get the mesh values
npts = 1000; % needs to be a power of 2
xi = linspace(xmin,xmax,npts);


%% bandwidth estimation

   % revised: sometimes KDE performs very poorly; a simple rule is just to
   % take the estimated bandwidth that is closest to range / 20
   bw_kde = kde(data_vec,npts,xmin,xmax);
   [xxx, xxx, bw_matlab] = ksdensity(data_vec,xi);
  
   %tar = range / 20;
   
   %bw_geo = geomean([bw_kde bw_matlab]);
   
   % RJE rule
   
   if bw_kde > 1 && bw_matlab > 1
       % take the smaller one
       bw_final = min(bw_kde, bw_matlab);
   else
       % bw_kde < 1 && bw_matlab < 1 or one is < 1 and the other is > 1
       % take the larger one
       bw_final = max(bw_kde, bw_matlab);
   end
   
% bw_all = nan(ncol,1);
% bw_vs_range = nan(ncol,1);
% 
% for m = 1:ncol
%    % revised: sometimes KDE performs very poorly; a simple rule is just to
%    % take the estimated bandwidth that is closest to 1.0
%    bw_kde = kde(data(:,m),npts,xmin,xmax); 
%    bw_matlab = ksdensity(data(:,m),xi);
%    
%    bw_vs_range(m) = (max(data(:,m)) - min(data(:,m))) / bw_all(m);
% end
% 
% %get the minimum observed bw value
% bw_final = min(bw_all)
% 
% % RJE manual adjustment
% rat = 50;
% if range / bw_final > rat % the original range of the full data set
%     bw_final = range / rat
% end

% updated method: just do kde on the *unique* values in the full dataset!
% no, this doesn't work well
%bw_final = kde(unique(data(:)),npts,xmin,xmax)

%% now use ksdensity.m 
f = nan(npts,ncol);


for m = 1:ncol
   fvals  = ksdensity(data(:,m),xi,'width',bw_final); % use the same bandwidth value for *all* columns
   f(:,m) = fvals(:);
   
   %if plot_figs == 1
       % plot result
       
       figure(200)
       subplot(2,1,1)
       hist(data_vec,20)
       set(gca,'xlim',[xmin xmax])
       
       subplot(2,1,2)
       plot(xi,fvals(:))
       hold on
   %end
   
end

hold off
set(gca,'xlim',[xmin xmax])

%% calculate the minimum overlap of ALL unique pairs

ovl_mat = nan(ncol,ncol);
KSstat = nan(ncol,ncol);
KSpval = nan(ncol,ncol);
ADstat = nan(ncol,ncol);
ADpval = nan(ncol,ncol);

% get the unique pairs
pairs = nchoosek(1:ncol,2);

npairs = size(pairs,1);

for p = 1:npairs
   f1 = f(:,pairs(p,1));
   f2 = f(:,pairs(p,2));
   
   ff = [f1 f2];
   fmin = min(ff,[],2); % f1 and f2 will have the same x-values
   
   if plot_figs == 1
       figure(201)
       plot(xi,ff)
       hold on
       plot(xi,fmin,'r')
       hold off
   end
%    trapz(xi,f1)
%    trapz(xi,f2)
   
   ovl_mat(pairs(p,2),pairs(p,1)) = trapz(xi,fmin); % the overlap coefficient is simply this area, which will have a range of 0 to 1
   
   %% other statistics
   [xxx KSp KSd] = kstest2(f1,f2); 
   
   KSstat(pairs(p,2),pairs(p,1)) = KSd;
   KSpval(pairs(p,2),pairs(p,1)) = KSp;
   
   inds = [ones(numel(f1),1); ones(numel(f2),1)+1];
   f3 = [f1; f2];
   X = [f3 inds];
   
   % AD test, adjusted for ties
   ADout = AnDarksamtest(X,.05);
   ADstat(pairs(p,2),pairs(p,1)) = ADout.ties_stat;
   ADpval(pairs(p,2),pairs(p,1)) = ADout.ties_pval;   
   
end


%% now we need to find the critical value using bootstrap analysis
% have to do this for each PAIR independently (permutation testing)

ovals = nan(iter,ncol); 

progressbar(0,0) % initialize
 
for m = 2:ncol
    progressbar((m-1)/ncol,0)
    for i = 1:iter
        progressbar([],i/iter)
        % pull a single value randomly from each *row*
        %c1 = randint(nsub,1,[1 ncol]);
        %c2 = randint(nsub,1,[1 ncol]);

        % get the real data        
        r = [data(:,1) data(:,m)];

        c = randint(nsub,1,[1 2]); % unique for each iter

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

        % ksdensity using the same bandwidth as used above
        f1 = ksdensity(p1,xi,'width',bw_final);
        f2 = ksdensity(p2,xi,'width',bw_final);

        f1 = f1(:);
        f2 = f2(:);

        % calculate the overlap
        ff = [f1 f2];
        fmin = min(ff,[],2); % f1 and f2 will have the same x-values

        if plot_figs == 1
           figure(202)
           plot(xi,ff)
           hold on
           plot(xi,fmin,'r')
           hold off
        end

        ovals(i,m) = trapz(xi,fmin); % so we have column 1 just be NaN

    end % iter loop
    
end
progressbar(1)

%% now KDE for the OVL values

xmax = 1.0; % easy
xmin = min(ovl_mat(:)); % from the above calcualtions
range = xmax - xmin;
%xmin = xmin - range/range_div;
xmin = 0.5;
xi_crit = linspace(xmin,xmax,npts);

% % get the minimum bandwidth as above
% bw_crit = nan(ncol,1);
% 
% for m = 2:ncol
%    bw_crit(m) = kde(ovals(:,m),npts,xmin,xmax); 
% end
% 
% bw_crit_min = min(bw_crit);

% just force this
bw_crit_min = .005;

% save the crit vals
x_crit_kde = nan(ncol,1);
x_crit_prc = nan(ncol,1);

for m = 2:ncol
    
    ovals_cur = ovals(:,m); 
    
    % now get the PDF    
    f_crit  = ksdensity(ovals_cur,xi_crit,'width',bw_crit_min);

    % now get the CDF
    cdf_crit  = ksdensity(ovals_cur,xi_crit,'width',bw_crit_min,'function','cdf');

    % and get the one-tailed critical value (p = .05)
    x_crit_prc(m) = prctile_nist(ovals_cur,5);
    x_crit_kde(m) = max(xi_crit(cdf_crit <= .05));

    figure(203)
    subplot(2,1,1)
    hist(ovals_cur,20)
    hmax = hist(ovals_cur,20);
    hmax = max(hmax);
    
    hold on
    plot([x_crit_prc(m) x_crit_prc(m)],[0 10],'r')
    set(gca,'xlim',[xmin xmax])
    hold off
    
    subplot(2,1,2)
    plot(xi_crit,f_crit,'LineWidth',2)
    hold on
    plot([x_crit_kde(m) x_crit_kde(m)],[0 max(f_crit)],'r')
    hold off
    pause(1)
end


%% output
output.data_bw_kde = bw_kde;
output.data_bw_matlab = bw_matlab;
output.data_bw_final = bw_final;
output.data_ovl = ovl_mat;

output.perm_bw = bw_crit_min;
output.perm_ovl = ovals;

output.x_crit_kde = x_crit_kde;
output.x_crit_prc = x_crit_prc;
output.ovl_gt = ovl_mat(:,1); % versus ground truth (column 1)
output.sig_kde = output.ovl_gt < x_crit_kde;
output.sig_prc = output.ovl_gt < x_crit_prc;

output.KSstat = KSstat(:,1);
output.KSpval = KSpval(:,1);
output.ADstat = ADstat(:,1);
output.ADpval = ADpval(:,1);



