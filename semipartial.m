function [output] = semipartial

% 1. Entr the data vector
data = input(' Enter an [N x 1] vector of dependent values, in [ ]: ');
nsub = size(data,1);

% 2. Enter the regressor matrix
regr = input(' Entr an [N x M] matrix of regressor values, in [ ]: ');
nreg = size(regr,2);

% check size
if size(data,1) == size(regr,1);
   % OK
   
else
    return
end

% for regression module to work in matlab, we must have a vector of 1s in
% the regressor matrix; i.e., the intercept
% see: http://www.mathworks.com/help/toolbox/stats/regress.html

regr(1:nsub,nreg+1) = ones(nsub,1); 

% 3. Do jackknife?
opt = input(' [1] Standard analysis; [2] delete-1 jackknife; [3] Monte Carlo permutations: ');

if opt == 1
   do_perm = 0;
   nperm = 1;
elseif opt == 2
   do_perm = 0;
   nperm = nsub;
elseif opt == 3
   do_perm = 1;
   nperm = input(' Number of permutations (default = 5000): ');
end

semipart_corrs = zeros(nperm,nreg); % individual runs as appropriate

% reformat the data

if opt == 1
   % just do the full partial correlation analysis 
   
   for m = 1:nreg
   y = 1:nreg;
   y = y(y~=m);
   
   X = [data regr(:,m)];
   Y = regr(:,y); % don't include the regressor just put into X
   
   pc = partialcorr(X,Y);
   part_corrs(1,m) = pc(1,2); % just the correlation value
           
   end
   
elseif opt == 3
   
   for p = 1:nperm
   
   % first, get a permutation of the data (not the regressors)
   
   xperm = randperm(nsub); % a random permutation of integers 1:nsub
   datap = data(xperm); % a reordering of the data
   
       for m = 1:nreg
       y = 1:nreg;
       y = y(y~=m);

       X = [datap regr(:,m)]; % now we use the permuted vector
       Y = regr(:,y); % don't include the regressor just put into X

       pc = partialcorr(X,Y);
       semipart_corrs(p,m) = pc(1,2); % just the correlation value

       end
   
   end
    
elseif opt == 2 

    for n = 1:nsub
        
    x = 1:nsub;
    x = x(x~=n); 
    
    data2 = data(x,:); % exclude one subject
    regr2 = regr(x,:);

       for m = 1:nreg
           y = 1:nreg;
           y = y(y~=m);

           X = [data2 regr2(:,m)];
           Y = regr2(:,y); % don't include the regressor just put into X

           pc = partialcorr(X,Y);
           semipart_corrs(n,m) = pc(1,2); % just the correlation value

       end
    end 
    
end

% convert to Fisher z-values
r = part_corrs;
z = 0.5 .* log((1+r)./(1-r));
       
% summary statistics
output.r_mean = mean(r);
output.r_std = std(r);
output.r_p95 = prctile(r,95);
output.r_p05 = prctile(r,5);

output.z_mean = mean(z);
output.z_std = std(z);
output.z_p95 = prctile(z,95);
output.z_p05 = prctile(z,5);

output.r_vals = r;
output.z_vals = z;


