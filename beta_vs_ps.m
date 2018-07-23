function out = beta_vs_ps(a1,b1,a2,b2,subN,iter,p_crit)

if nargin < 5
    subN = 20:20:100;
    
else
    if isempty(subN)
        subN = 20:20:100;
    end
        
end

if nargin < 6
    iter = 20000;
end

if nargin < 7
    p_crit = .05;   
end

% how many values do we need?
minV = max(subN) * iter;

if minV < 1000000
    minV = 1000000;
end


% simulate full set of values
plot_it = 1;
sim = beta_2sim(a1,b1,a2,b2,minV,plot_it);

% preallocate
ps_mat = nan(iter,numel(subN));
p_mat  = nan(iter,numel(subN)); % p-value from MW test

progressbar(0,0)

for j = 1:numel(subN)
    
    % index locations - do up front to save time
    ind_all = randperm(minV);
    
    % reshape the index matrix
    nrow = subN(j);
    ncol = iter;
    
    ind = reshape(ind_all(1:nrow*ncol),nrow,ncol);
    
    progressbar((j-.01)/numel(subN),[])
       
    for i = 1:iter
        
        % get the values
        sub1 = sim.sim1(ind(:,i));
        sub2 = sim.sim2(ind(:,i));
        
        % PS statistic (samp2 vs samp1)
        rs = ranksum_rje(sub1,sub2);
        ps_mat(i,j) = rs.PS;
        p_mat(i,j)  = rs.Pval_approx;
        
        if rem(i,100) == 0
            progressbar([],i/iter)
        end

    end

end

progressbar(1)

% calculate small, medium, large effects
% Thresholds per Vargha 2000

abs_ps_mat = abs(ps_mat - 0.5);

% only count scores where p < .05!
p_inc = p_mat < p_crit; % logical

p_sig = 100 * mean(p_inc)

check_these = abs_ps_mat .* p_inc;

prc_small  = 100 * mean(check_these >= 0.06)
prc_med    = 100 * mean(check_these >= 0.14)
prc_large  = 100 * mean(check_these >= 0.21)

prc_gt10   = 100 * mean(check_these >= 0.10)

%% figures
fignum = 410;
boxplot_rje(ps_mat,[1 99],fignum);

% rescale y
ymin = min(prctile(ps_mat,1));
ymax = max(prctile(ps_mat,99));
ylim([ymin ymax])

figure(420)
for j = 1:numel(subN)
    plot(ps_mat(:,j),log10(p_mat(:,j)),'.','DisplayName',['N = ' num2str(subN(j))])
    hold on
end
hold off
xlabel('PS statistic')
ylabel('log10(P)')
legend('show','Location','northeastoutside')
xlim([0 1])

%% outputs
out.minV = minV;
out.true_pr_2gt1 = sim.pr_2gt1;
out.subN         = subN;
out.ps_mat       = ps_mat;
out.p_mat        = p_mat;
