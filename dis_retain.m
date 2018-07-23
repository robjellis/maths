% [retained] = dis_retain(X, d)
%
% evaluate all pairs of XYZ coordinates (where X is an N-by-3 matrix
% and provide an N-by-1 list (1 vs 0) of
% coordinates that are at least d units apart.
%
% rje 2012

function [retained] = dis_retain(X,d)

% how big?

if size(X,2) ~= 3
    % problem
    fprintf('\n Error: only an N-by-3 input matrix is accepted.\n')
    return
end

N = size(X,1); % number of coordinates

temp = zeros(N,N-1);
chk = 1:N;

for i = 1:N
   
   a = X(i,:);

       chkr = chk(chk~=i); % the remaining
       for j = 1:numel(chkr)
           b = X(chkr(j),:);

           % now we do the calculation
           x1 = a(1); y1 = a(2); z1 = a(3);
           x2 = b(1); y2 = b(2); z2 = b(3);

           dis = sqrt((x2 - x1)^2 + (y2 - y1)^2 + (z2 -z1)^2);
           
           % is it greater than min?
           temp(i,j) = dis >= d;
       end
   
end

temp

% now just sum the columns

temp2 = sum(temp,2);

retained = temp2 == N-1
