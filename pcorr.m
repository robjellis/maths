function [vals partialr] = pcorr(y,X,var)

% partial correlation for variable VAR
% [vals partialr] = pcorr(y,X,var)

% need to eliminate any data that have NaN


   xx = 1:size(y,1); % index
   xx = xx(:);
   %size(isnan(data))
   yy = xx(isnan(y)==0); % retain these rows for data and regr
   IDs = (1:size(y,1))'; % case numbers for reference
   
   y = y(yy);
   X = X(yy,:); % take from the ORIGINAL regressor matrix
   IDs = IDs(yy); % retained cases
   
   nsub = size(y,1); % redefine this now
   
   nreg = size(X,2);
   
% ***** for regression module to work in matlab, we must have a vector of 1s in
% the regressor matrix; i.e., the intercept
% see: http://www.mathworks.com/help/toolbox/stats/regress.html

Xmod = X;
Xmod(1:nsub,nreg+1) = ones(nsub,1); % add it at the end
k = 1:nreg+1;

% pull out the variable we want
   
   x = Xmod(:,var);
   keep = k(k~=var);
   Xmod = Xmod(:,keep);
   
   % now we get the residuals from (1) y on Xmod and (2) x on Xmod
   
   [b,bint,residx] = regress(x,Xmod);
   [b,bint,residy] = regress(y,Xmod);
   
   %figure
   %plot(residx,residy,'.')
   partialr = corr(residx,residy); % rje confirmed this is correct
   vals = [residx residy IDs];
   
   
   
