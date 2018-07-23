function [r rr] = rje_pearson(x,y,p,figs)

% rje_pearson(x,y,p,figs) 
% pearson correlation that does not require matlab stats toolbox
%
% operates on a single [N x 1] X vector and [N x 1] Y vector
%
% rje, august 2012

% note: MATLAB corrcoef.m will do the r-value without use of the stats toolbox,
% but won't do the "robust" operation performed here.

% check size
x = x(:);
y = y(:);
sx = size(x);
sy = size(y);

n = sx(1);

if sum(sx == sy)==2 && sx(2)==1

  % http://en.wikipedia.org/wiki/Pearson_product-moment_correlation_coefficient

  cv = cov(x,y); % note: the default will use the N - 1 correction for unbiased estimation of a sample
  cv = cv(1,2);
  r = cv / (std(x)*std(y)); 
  
  % r = (1 / (n - 1)) * sum( (x-mean(x))/std(x) .* (y-mean(y))/std(y) );
  
  % out of curiosity, what happens if we trim the tails and then do the
  % correlation?
  

  
  % special case for p = 0 ... include all values
  if p == 0
      xx = x;
      yy = y;
  else
    prs = [p 100-p];
    px = rje_spmprc(x,prs);
    py = rje_spmprc(y,prs);
        
    xx = x.*(x>px(1)); xx = xx.*(xx<px(2));
    yy = y.*(y>py(1)); yy = yy.*(yy<py(2));
  end
  zz = (xx~=0) + (yy~=0); % only those pairs that are non-zero
  
  xx = xx(zz==2);
  yy = yy(zz==2);
  nn = size(xx,1);
  
  cv2 = cov(xx,yy);
  cv2 = cv2(1,2);
  rr = cv2 / (std(xx)*std(yy));
  %rr = (1 / (nn - 1)) * sum( (xx-mean(xx))/std(xx) .* (yy-mean(yy))/std(yy) )
  
  % figures
  if figs == 1
     figure(100)
     plot(x,y,'.')
     hold on
     plot(xx,yy,'r.')
     hold off
  end
else
  fprintf(' Warning: "pearson" operates on a single X and Y vector only. Respecify.\n');
  r = [];
  rr = [];
  return
end

