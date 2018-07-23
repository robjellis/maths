function sim_res = sampling_dist(X,stat,samp,iter,saveit)

% get the sempling distribution of STAT statistic from X distribution
% stat can be 
%             - 'p50': the median
%             - 'p75': 75th percentile
%             - 'tm': trimean
%             - 'thr': custom threshold derived from the distribution
%

if isempty(saveit)
    saveit = input(' Save result? \n   [1] No \n   [2] Yes \n   --> ');   
end

if saveit == 1
    % do nothing
elseif saveit == 2
    name = input(' Name of *.mat file: ','s');
end

tic;

X = X(:);
nX = numel(X);

% do before the loop
if strcmp(stat,'thr')
    thr = prctile(X,50);
    prop_above_thr = nan(iter,1);
else
    thr = [];
%    rep_rate_med = nan(iter,1);
%    rep_rate_max = nan(iter,1);
end


%%%%%%%%
% must do iteratively for sake of memory

stat_res = nan(iter,1);

progressbar(0)
for i = 1:iter
    inds = randint(samp,1,[1 nX]);
    vals = X(inds);

    if strcmp(stat,'thr')
        vals_keep = vals; % make a copy
        vals_keep(vals <= thr) = NaN; % turn into NaN so prctile not affected

        vals_count = vals > thr; % 1 or 0
        prop_above_thr(i) = sum(vals_count) ./ samp; % proportion in each row     
        
        stat_res(i) = prctile(vals_keep,50);
    elseif strcmp(stat,'p50')
        stat_res(i) = prctile(vals,50);
    elseif strcmp(stat,'p75')
        stat_res(i) = prctile(vals,75);
    elseif strcmp(stat,'tm')
        stat_res(i) = trimean(vals);
        
%         % other stats
%         T = tabulate(vals(i,:));
%         rep_rate_med(i) = median(T(:,3));
%         rep_rate_max(i) = max(T(:,3));
    end
    
    progressbar(i/iter)
end

progressbar(1)


duration = toc;

figure(100)
ecdf(stat_res)
xlabel('ICF percentile')


%% start writing results
sim_res.stat = stat;
sim_res.iter = iter;
sim_res.samp = samp; 
sim_res.thr  = thr;
sim_res.duration = duration;
sim_res.min = min(stat_res);
sim_res.md  = median(stat_res);
sim_res.max = max(stat_res);
sim_res.skewness = skewness(stat_res,0);
sim_res.stat_res = stat_res;

%% figures
if strcmp(stat,'thr')
    
   figure(103)
   ecdf(prop_above_thr)
   
   sim_res.prop_above_lo = prctile(prop_above_thr,2.5);
   sim_res.prop_above_hi = prctile(prop_above_thr,97.5);   
   
else
%     figure(101)
%     subplot(1,2,1)
%     ecdf(rep_rate_med)
%     xlabel('Repetition rate median (%)')
%     ylabel('Cumulative prob.')
% 
%     subplot(1,2,2)
%     ecdf(rep_rate_max)
%     xlabel('Repetition rate maximum (%)')
%     ylabel('Cumulative prob.')    
end

% figure(101)
% qqplot(res)



if saveit == 2
   expr = ['save ' name '.mat sim_res'];
   eval(expr)
end

