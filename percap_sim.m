function [t_test] = percap_sim(nsub,nvox,nsim,upper)

% 
%function [t_test] = percit_sim(nsub,nvox,nsim,upper)
%
% simulation to see whether we get Student T distributions data with small N percap
% transformations
%
% upper ranges from 0 to 100 (e.g., 99) is the 99th percentile)


% mean adjustment?

meanadj = input(' Apply a mean parameter to randn?: \n  [1] all means = 0; \n  [2] fixed mean for all subjects; \n  [3] random mean (between -0.5 and 0.5): ');

% get lower
lower = (100 - upper);

% number of subjects

rdata = zeros(nsub,nvox);
pdata = zeros(nsub,nvox);

pvals = [.05 .01 .005 .001 .0005];
lpvals = -1*log10(pvals);
prvals = [95 99 99.5 99.9 99.95]; 

nump = numel(pvals);

pc_sim = zeros(nsim,nump);
rc_sim = zeros(nsim,nump);
diff_sim = zeros(nsim,nump);

% counter
ctr = floor(nsim * [.1 .2 .3 .4 .5 .6 .7 .8 .9 1.0]);

k = 1;

if meanadj == 2
   madj2 = input(' Enter a fixed mean for randn (e.g., 0.3): '); 
end

fprintf('\n Working ');


    
for j = 1:nsim
    
    if j == ctr(k)
       fprintf('.')
       if k < 10
            k = k + 1;
       end
    end
    
    for i = 1:nsub
        
        x = randn(1,nvox);
        
        meanx = mean(x);
        
        x = x - meanx; 
        
        if meanadj == 1
            x = x - meanx; % to make sure distribution is at 0
        elseif meanadj == 2
            x = x - meanx + madj2;
        elseif meanadj == 3
            x = x - meanx + rand - 0.5; % random number between -0.5 and 0.5 
        end
        
        
        
        [percvol] = percap(x,upper);

        pdata(i,:) = percvol;
        rdata(i,:) = x;

    end

    % do the t-test on percap
    [a b c d] = ttest(pdata);
    p_tvals = d.tstat; % these are correct t-values; checked against Statistica (don't care 1- or 2-tailed p-values right now)

    % do it on the original data
    [a b c d] = ttest(rdata);
    r_tvals = d.tstat;

    % histogram of the percap t-distribution
    %figure(50)
    %hist(p_tvals,50);


    df = nsub - 1;

    % theoretical values
    tcrit = abs(tinv(pvals,df));

    % observed critical values from random data
    rcrit = prctile(r_tvals,prvals);

    % observed critical values from percap data
    pcrit = prctile(p_tvals,prvals);

    % collect for the M.C. simulation
    rc_sim(j,:) = rcrit;
    pc_sim(j,:) = pcrit;
    diff_sim(j,:) = pcrit - rcrit;
    
    if nsim == 1   
        figure(60)
        plot(tcrit,tcrit,'LineWidth',2,'Color',[0 0.498 0],'Marker','o');
        hold on
        plot(tcrit,rcrit,'bo-')
        plot(tcrit,pcrit,'ro-')
        hold off
        xlabel('theoretical T-crit')
        ylabel('(observed) T-crit')
    end
        
end

% diff
figure(61)
hist(diff_sim,50);

% paired t-test
[a b c d] = ttest(diff_sim);

t_test.p_crit = pvals;
t_test.t_val = d.tstat;
t_test.p_val = b;

% M.C. results

rc_mn = mean(rc_sim);
rc_sdp = rc_mn + std(rc_sim);
rc_sdm = rc_mn - std(rc_sim);

pc_mn = mean(pc_sim);
pc_sdp = pc_mn + std(pc_sim);
pc_sdm = pc_mn - std(pc_sim);

figure(70)

plot(tcrit,tcrit,'LineWidth',3,'Color',[0 0.498 0],'Marker','o');

hold on
plot(tcrit,rc_mn,'bo-')
plot(tcrit,rc_sdp,'Marker','o','LineWidth',2,'Color','b','LineStyle','--')
plot(tcrit,rc_sdm,'Marker','o','LineWidth',2,'Color','b','LineStyle','--')

plot(tcrit,pc_mn,'ro-')
plot(tcrit,pc_sdp,'r.','LineStyle','--')
plot(tcrit,pc_sdm,'r.','LineStyle','--')



        xlabel('theoretical T-crit')
        ylabel('(observed) T-crit')
        title(['N = ' num2str(nsub) '; L = ' num2str(lower) '; U = ' num2str(upper) '; ' num2str(nsim) ' Monte Carlo iterations']);
hold off


fprintf('\n\n Finished. \n\n');