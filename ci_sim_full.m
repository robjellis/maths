function out = ci_sim_full(N,width,iter)

%
% ci_sim_full(N,width,iter)
%
% use ci_sim.m to test the full space of p (success rates) and n sizes
%
% RJE | 2017.11.18

% clear figures
figure(100)
clf

figure(101);
clf

pause(0.2)

% defaults
P = 0:.005:1; % all methods tested have a manual correction for x == 0 or x == n
nmeas = 4; 

N_orig = N; % keep this

if N == 0 % just a simple way of doing this
    N = [8 16 32 64 128 256 512 1024 2048];
end

if nargin < 2
    width = 95;
end

if nargin < 3
    iter = 1000;
end

numN = numel(N);
numP = numel(P);

% plot all or just some?
if max(N_orig) == 0 % just plot a key subset
    plot_N = [8 64 512]; % 2^3, 2^6, 2^9
else
    plot_N = N_orig; % the full set
end

% variables - currently at 4
cov_wald_nc = nan(numP,numN);
cov_wald_cc = nan(numP,numN);
cov_wils_nc = nan(numP,numN);
cov_wils_cc = nan(numP,numN);

% summary percentile stats (3 percentile values)
stats_wald_nc = nan(numN,3);
stats_wald_cc = nan(numN,3);
stats_wils_nc = nan(numN,3);
stats_wils_cc = nan(numN,3);

% start progress bar
progressbar(0,0)

% THE LOOP
for n = 1:numN
    progressbar(n/numN - eps,[])
    
    for p = 1:numP
        progressbar([],p/numP - eps)
    
        res = ci_sim(P(p),N(n),width,iter);
        
        % store the variables
        cov_wald_nc(p,n) = res.inc_wald_nc;
        cov_wald_cc(p,n) = res.inc_wald_cc;
        cov_wils_nc(p,n) = res.inc_wils_nc;
        cov_wils_cc(p,n) = res.inc_wils_cc;
        
    end
    
    
    % *********************************************** 
    % Range from low to hi percentile
    % set the width == the main parameter width
    
    lo = (100 - width)/2;
    hi = 100 - lo;
    med = 50;

    % original version
    stats_wald_nc(n,:) = prctile_nist(cov_wald_nc(:,n),[lo med hi])'; % transpose into a row vector
    stats_wald_cc(n,:) = prctile_nist(cov_wald_cc(:,n),[lo med hi])';
    stats_wils_nc(n,:) = prctile_nist(cov_wils_nc(:,n),[lo med hi])';
    stats_wils_cc(n,:) = prctile_nist(cov_wils_cc(:,n),[lo med hi])';
    
%     updated: do 100 - value so we can have better visual scaling; i.e., we are
%     plotting the % of cases in which the CI does *NOT* cover the true value
%     
%     stats_wald_nc(n,:) = prctile_nist(100 - cov_wald_nc(:,n),[lo med hi])'; % transpose into a row vector
%     stats_wald_cc(n,:) = prctile_nist(100 - cov_wald_cc(:,n),[lo med hi])';
%     stats_wils_nc(n,:) = prctile_nist(100 - cov_wils_nc(:,n),[lo med hi])';
%     stats_wils_cc(n,:) = prctile_nist(100 - cov_wils_cc(:,n),[lo med hi])';    
    
    % *********************************************** 
    % Percentage coverage figure
    if ismember(N(n),plot_N)
        for m = 1:nmeas
            if m == 1
                var = cov_wald_nc;
                tit = 'Wald, no correc.';
            elseif m == 2
                var = cov_wald_cc;
                tit = 'Wald, cont. correc.';
            elseif m == 3
                var = cov_wils_nc;
                tit = 'Wilson, no correc.';
            elseif m == 4
                var = cov_wils_cc;
                tit = 'Wilson, cont. correc.';
            end

            % ****************
            figure(100)
            subplot(2,2,m)
            plot(P,var(:,n),'DisplayName',['n = ' num2str(N(n))])
            hold on

            title(tit)
            xlabel('True P')
            ylabel('% coverage of true P')

            % set the axes
            ylim([width-15, 100]) % or 100.1 so we can see the top of the error bar
            set(gca,'box','off')

            if m == 4 && N(n) == max(plot_N) % wait until the end otherwise it won't show all values
               % note: there is no function to add a "global" legend when using subplot
               legend('show','Location','southeast');
            end

        end % m loop
    end
end % n loop

progressbar(1)

% *********************************************** 
% Coverage consistency figure (figure101)
figure(101)    

for m = 1:nmeas
    if m == 1
        var = stats_wald_nc;
        tit = 'Wald, no correc.';
    elseif m == 2
        var = stats_wald_cc;
        tit = 'Wald, cont. correc.';
    elseif m == 3
        var = stats_wils_nc;
        tit = 'Wilson, no correc.';
    elseif m == 4
        var = stats_wils_cc;
        tit = 'Wilson, cont. correc.';
        
    end

    % ****************
    % here we carefully set up the x-axis so we have nice log scaling!
    subplot(2,2,m,'XMinorTick','on','XScale','log','XTick',N); % same values as input N
    hold('on') % keep this here
    
    % error bars are *relative* to the median
    med = var(:,2);
    neg  = med - var(:,1);
    pos  = var(:,3) - med;
    
    % https://www.mathworks.com/help/matlab/ref/errorbar.html
    % note: the end caps are small, but there is no way to edit this simply.
    errorbar(N,med,neg,pos,'-o','MarkerSize',5,'MarkerEdgeColor',[0 0.45 0.75],'MarkerFaceColor',[0 0.45 0.75])     

    % set the axes
    xlim([min(N)/2 max(N)*2])
    ylim([width-15, 100])
    title(tit)
    xlabel('n')
    ylabel('% coverage of true P')
    set(gca,'box','off') % makes it easier to see the upper bound of the C.I.

end % m loop

% store all the stats
stats.wald_nc = stats_wald_nc;
stats.wald_cc = stats_wald_cc;
stats.wils_nc = stats_wils_nc;
stats.wils_cc = stats_wils_cc;

hold off








%% outputs
out.stats = stats;