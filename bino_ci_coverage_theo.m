function out = bino_ci_coverage_theo(width,N)

if nargin < 1
    width = 95;
end

if nargin < 2
    N = [8 16 32 64 128 256 512 1024];
end

% turn off warnings for large binomial coefficents
warning('off','MATLAB:nchoosek:LargeCoefficient');

beep off

% test these
P = 0:.005:1;

numN = numel(N);
numP = numel(P);

% main storage
cov_wald_nc = nan(numP,numN);
cov_wald_cc = nan(numP,numN); 
cov_wils_nc = nan(numP,numN); 
cov_wils_cc = nan(numP,numN); 

progressbar(0,0)
fprintf([' Computing theoretical coverage for ' num2str(numel(N)) ' sample size values ... ']);

for n = 1:numN
    progressbar(n/numN-.01,[])
    
    numX = N(n)+1; % to count 0
    
    for p = 1:numP

        progressbar([],p/numP)
        
        % temporary storage - reset each time
        C = nan(numX,1);
        
        ind_wald_nc = zeros(numX,1);
        ind_wald_cc = zeros(numX,1);
        ind_wils_nc = zeros(numX,1);
        ind_wils_cc = zeros(numX,1);
        
        
        for x = 0:N(n) % store values for each value of x

            % from the equation        
            C(x+1,1) = nchoosek(N(n),x) * P(p)^x * (1 - P(p))^(N(n)-x); % store it            

            % the CIs for each method
            B = bino_ci_calc(x,N(n),width);

            if B.wald_nc(1) <= P(p) && P(p) <= B.wald_nc(2) 
                ind_wald_nc(x+1) = 1;
            end

            if B.wald_cc(1) <= P(p) && P(p) <= B.wald_cc(2)
                ind_wald_cc(x+1) = 1;
            end

            if B.wils_nc(1) <= P(p) && P(p) <= B.wils_nc(2)
                ind_wils_nc(x+1) = 1;
            end   

            if  B.wils_cc(1) <= P(p) && P(p) <= B.wils_cc(2)
                ind_wils_cc(x+1) = 1;
            end
            
        end % x loop
        
        % get the sum of the product
        cov_wald_nc(p,n) = sum(C .* ind_wald_nc);
        cov_wald_cc(p,n) = sum(C .* ind_wald_cc);
        cov_wils_nc(p,n) = sum(C .* ind_wils_nc);
        cov_wils_cc(p,n) = sum(C .* ind_wils_cc);
        
    end % p loop
    
end % N loop
        
progressbar(1) 
fprintf(' Done.')


figure(30)
clf

figure(31)
clf

for m = 1:4
    
    if m == 1
        data = cov_wald_nc;
        tit = 'Wald, no correc.';
    elseif m == 2
        data = cov_wald_cc;
        tit = 'Wald, cont. correc.';
    elseif m == 3
        data = cov_wils_nc;
        tit = 'Wilson, no correc.';
    elseif m == 4
        data = cov_wils_cc;
        tit = 'Wilson, cont. correc.';
    end
    
    
    % for figure 31
    med = nan(1,numel(N));
    lb  = nan(1,numel(N)); % lower bound of interval
    ub  = nan(1,numel(N)); % upper bound of interval
    
    for n = 1:numN
        
        this_data = data(:,n);
        
        % plot coverage as a function of P and N
        figure(30)
        subplot(2,2,m)
        plot(P,this_data)
        title(tit)
        ylabel('Coverage of true P')
        xlabel('P')
        hold on
        
        % summary statistics
        med(n) = median(this_data);
        prc_lo = (100-width)/2;
        prc_hi = 100 - prc_lo;
        lb(n)  = prctile_nist(this_data,prc_lo); % lower bound of interval
        ub(n)  = prctile_nist(this_data,prc_hi);    % upper bound of interval
    
    end
    
    figure(30)
    ylim([width/100 - .15 1.0])
    
    % **************** figure 31 ****************
    figure(31) % no hold on here ...
    
    % error bars are *relative* to the median
    neg  = med - lb;
    pos  = ub - med;

    % here we carefully set up the x-axis so we have nice log scaling!
    subplot(2,2,m,'XMinorTick','on','XScale','log','XTick',N); % same values as input N
    hold('on') % keep this here

    % https://www.mathworks.com/help/matlab/ref/errorbar.html
    % note: the end caps are small, but there is no way to edit this simply.
    errorbar(N,med,neg,pos,'-o','MarkerSize',5,'MarkerEdgeColor',[0 0.45 0.75],'MarkerFaceColor',[0 0.45 0.75])     

    % set the axes
    xlim([min(N)/2 max(N)*2])
    ylim([width/100 - .15 1.0])
    title(tit)
    xlabel('n')
    ylabel('Coverage consist. across P')
    set(gca,'box','off') % makes it easier to see the upper bound of the C.I.
   
end % subplot (m) loop

figure(30)
hold off

figure(31)
hold off


% ****** histograms - too messy
%bins = 30;
% figure(35)
% subplot(2,2,1)
% hist(cov_wald_nc,bins)
% xlim([.8 1.0])
% 
% subplot(2,2,2)
% hist(cov_wald_cc,bins)
% xlim([.8 1.0])
% 
% subplot(2,2,3)
% hist(cov_wils_nc,bins)
% xlim([.8 1.0])
% 
% subplot(2,2,4)
% hist(cov_wils_cc,bins)
% xlim([.8 1.0])

out = [];
