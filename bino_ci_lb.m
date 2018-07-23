function bino_ci_lb(width,N)

%
% bino_ci_lb(width,N)
%
% Compute the lower bound of binomial CIs of four methods
% width: e.g., 80, 90, 95
%
% RJE | 2017.11.18

if nargin == 0
    width = 95;
end

if nargin < 2
    N = [8 16 32 64 128 512 1024];
end

%P = 0:.01:1.0; % target values; will be converted into actual integers

numN = numel(N);

nmeas = 4;

figure(20) % leave this here
clf

figure(21)
clf

fprintf('\n Working ... ')
progressbar(0)
for n = 1:numN
    
    % get the success cases
    %X = unique(round(N(n).*P))
    X = 0:N(n); % do it manually, it's still fast; this way we evalulate every possible value!

    % get the observed proportion
    p = X ./ N(n);
    p = p';
    
    % how many?
    numX = numel(X);
    
    % store lower bounds
    LB = nan(numX,4); % four measures
    
    for x = 1:numX
        
        % get the CI lower bounds
        R = bino_ci_calc(X(x),N(n),width);
        
        LB(x,1) = R.wald_nc(1);
        LB(x,2) = R.wald_cc(1);
        LB(x,3) = R.wils_nc(1);
        LB(x,4) = R.wils_cc(1);
        
    end
    
    for m = 1:nmeas
        if m == 1
            tit = 'Wald, no correc.';
        elseif m == 2
            tit = 'Wald, cont. correc.';
        elseif m == 3
            tit = 'Wilson, no correc.';
        elseif m == 4
            tit = 'Wilson, cont. correc.';
        end

        figure(20)
        subplot(2,2,m)
        plot(p,LB(:,m),'DisplayName',num2str(N(n)))
        hold on
        
        if m == 1 && n == numN
           % add the legend; there is no function to add a "global" legend when using subplot
           legend('show','Location','northwest');
        end

        title(tit)
        xlabel('Observed p')
        ylabel('Lower bound of CI')
        ylim([0 1]) % consistent for all subplots
        
        % ***** now LB divided by p
        figure(21)
        subplot(2,2,m)
        plot(p,LB(:,m)./p,'DisplayName',num2str(N(n)))
        hold on
        
        title(tit)
        xlabel('Observed p')
        ylabel('LB / Observed p')
        ylim([0 1]) % consistent for all subplots
        
        if m == 4 && n == numN
           % add the legend; there is no function to add a "global" legend when using subplot
           legend('show','Location','southeast');
        end        
    end
    progressbar(n/numN)
end

progressbar(1)

fprintf(' Done.\n')

figure(20)
hold off

figure(21)
hold off

