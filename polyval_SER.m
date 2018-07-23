function [SER order] = polyval_SER(x,y,ord)

% standard error of regression for polynomials

x = x(:);
y = y(:);

if isempty(ord)
    ord_check = [1 2 3 4 5];
else
    ord_check = ord; % could be multiple values too
end

if nargin < 4
    plot_it = 1;
end

xmin = min(x);
xmax = max(x);
ymin = min(y);
ymax = max(y);

numx = numel(x);

% unique x-values for plotting
x_u = unique(x);

figure(25)
clf
plot(x,y,'.')
hold on
axis([xmin xmax ymin ymax])
   
SS_resid = nan(numel(ord_check),1);
SER      = nan(numel(ord_check),1);

progressbar(0,0)

for j = 1:numel(ord_check)
  progressbar(j/numel(ord_check),[])
  % check the model on the full data
  this_ord = ord_check(j);
  coef = polyfit(x, y, this_ord);

  % get all y-values for residual analysis
  y_fit = polyval(coef,x);
  
  % get y-values for plotting
  y_plot = polyval(coef,x_u);

  % plot it
  if plot_it == 1
      figure(25)
      plot(x_u,y_plot,'r','LineWidth',1.25)
      axis([xmin xmax ymin ymax])
  end
  
  yresid = y - y_fit;
  
  % the statistic we care aboutg is the sum of squares of residuals
  SS_resid(j) = sum(yresid.^2);
  
  % now calculate the SER
  SER(j) = sqrt(SS_resid(j) / (numx - ord_check(j) - 1));
  
end % j loop

progressbar(1)
  
if plot_it == 1
    figure(25)
    hold off
    
    figure(26)
    semilogy(ord_check,SS_resid,'Marker','.')
    xlabel('Polynomial order')
    ylabel('SS residual')
    
    figure(27)
    semilogy(ord_check,SER,'Marker','.')
    xlabel('Polynomial order')
    ylabel('Std Error of Regression (SER)')   
end

order = ord_check;

