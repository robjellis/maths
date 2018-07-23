function out = robustfit_vs_ols(x,y)

x = x(:);
y = y(:);

% compre OLS versus IRLS (robustfit.m)

% OLS
% add a vector of 1s to X
col = ones(numel(x),1);
X = [x col];

[ols_b ols_bint] = regress(y,X)

% IRLS
[irls_b, irls_stats] = robustfit(x,y)

%% fit model
xfit    = [min(x) max(x)];
y_ols   = ols_b(1)*xfit + ols_b(2);
y_irls  = irls_b(2)*xfit + irls_b(1);

%% PLOTS
figure(100)
plot(x,y,'.')
hold on
    plot(xfit,y_ols,'r','DisplayName','OLS')
    plot(xfit,y_irls,'b','DisplayName','IRLS')

hold off
xlabel('X values')
ylabel('Y values')
legend('show','Location','NorthEastOutside')