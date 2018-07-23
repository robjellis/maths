function [ranks] = rankit(x,rnd)

% retuns the ranks of a single [N x 1] vector. 
% Ties are appropriately handeled
% rnd = number of decimals to round to (0 to 4)
%
% by RJE

x = x(:); % just to be sure

% round it

if rnd == 0
   mult = 1;
elseif rnd == 1
   mult = 10;
elseif rnd == 2
   mult = 100;
elseif rnd == 3
   mult = 1000;
elseif rnd == 4
   mult = 10000;
end

x = round(x*mult)/mult;

numx = numel(x);

pos = (1:numx)'; % position
ranks = (1:numx)'; % actual ranks

vals = [pos x];

vals = sortrows(vals,2); % sort by x

tie = zeros(numx,1);

for i = 1:numx
    if i == 1
       if x(i + 1) == x(i)
          % stay at zero
       else
           tie(i) = 1;
       end
   elseif i < numx 
        if x(i - 1) == x(i) && x(i + 1) == x(i)
           % stay at zero
        elseif x(i -1) == x(i) && x(i + 1) > x(i)
           % stay at zero
        elseif x(i - 1) < x(i) && x(i + 1) == x(i)
           % stay at zero
        elseif x(i - 1) < x(i) && x(i + 1) > x(i)
            tie(i) = 1;
        end
    elseif i == numx
        if x(i - 1) == x(i)
            tie(i) = 0;
        else
            % stay at zero
        end
    end
end

tie
    

% are there ties?
diffx = diff(vals(:,2));
diffx = diffx ~=0; % non-ties


vals2 = [ranks vals]; 

vals3 = sortrows(vals2,2); % sort by pos to restore the original vector order

ranks = vals3(:,1);

