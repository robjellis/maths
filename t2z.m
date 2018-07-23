function [zvals] = t2z(tvals,df)

% [zvals] = t2z(tvals,df)
%
% takes a matrix of t-values and converts them to z-values using MATLAB
% functions
% 
% by rje

% update 2015.04.07; because matlab has trouble with numbers close to 1.0,
% but can more easily represent very small numbers (using exponents), we
% do the conversion by first taking the negative t-value and then restoring
% the sign

signs = sign(tvals);

tvals = abs(tvals);
tvals = -1* tvals; % this will always make p-values very small

pvals = tcdf_r14(tvals,df); % note: this is a direct match up of the two CDF functions
% this is a copy of the tcdf function from version 14 of matlab. Older
% versions have "missing" p-values (i.e., p = 0) for intermediate t-values;
% a strange situation which has been fixed in later versions of the script.

% get z
zvals = norminv(pvals);

% % Note: for some high d.f. values, Matlab doesn't calculate things correctly
% if min(zvals) == -Inf
%     % then we have missing values
%     
%     % get a full range
%     t = 0: -1 : - 50;
%     
%     p = tcdf(t,df);
%     
%     z = norminv(p);
%     
%     % use pchip interpolation
%     tt = 0 : - .0001 : -50;
%     
%     zz = pchip(t,z,tt);
%     
%     % find the locations
%     inds = find(zvals == -Inf);
%     
%     for i = 1:numel(inds)
%         % get the t-value
%         this_t = tvals(inds(i));
%         
%         % match it with z-value
%         this_z = max(zz(this_t > tt))
%         
%         % put it back
%         zvals(inds(i)) = this_z
%     end
% else
%     % everything is fine
% end


% now restore the correct sign
zvals = abs(zvals) .* signs;

