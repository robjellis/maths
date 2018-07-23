function out = beta_2sim(a1,b1,a2,b2,N,plot_it)

% simple way of comparing two beta distributions and get Pr(2>=1)
%
% RJE | 6 March 2018

if nargin < 5
    N = 1000000;
end

if nargin < 6
    plot_it = 1;
end

tic;

% theoretical PDFs
x = 0:.001:1;

pdf1 = betapdf(x,a1,b1);
pdf2 = betapdf(x,a2,b2);

% simulate values
sim1 = betarnd(a1,b1,N,1);
sim2 = betarnd(a2,b2,N,1);

% Pr(2>1)
pr_2gt1 = mean(sim2>=sim1); % 0 to 1 scale

tocc = toc;

%% plots

if plot_it == 1
    figure(400)
    clf

    subplot(1,2,1)
    plot(x,pdf1,'r','DisplayName','PDF for 1','LineWidth',1.2)  
    hold on
    plot(x,pdf2,'b','DisplayName','PDF for 2','LineWidth',1.0)  
    hold off    
    xlabel('Theoretical values')
    ylabel('Probability density')
    legend('show','Location','best')
    

    subplot(1,2,2)
    [f1, x1] = ecdf(sim1);
    [f2, x2] = ecdf(sim2);
    plot(x1,f1,'r','DisplayName','ECDF for 1','LineWidth',1.2) 
    hold on
    plot(x2,f2,'b','DisplayName','ECDF for 2','LineWidth',1.0) 
    hold off
    xlabel('Observed values')
    ylabel('Cumulative prob.')
    legend('show','Location','northwest')
end

set(gcf,'color','w');

out.mn_sim1 = mean(sim1);
out.mn_sim2 = mean(sim2);
out.pr_2gt1 = pr_2gt1;
%out.sim1 = sim1;
%out.sim2 = sim2;
out.duration_sec = tocc;