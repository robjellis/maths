function C = time_vec_align(A,B)

% take two vectors of timestamps A and B, and return adusted matrix C with the timestamps aligned as much as possible
%
% ** this code is not working correctly yet!

A = A(:);
B = B(:);

% size
numA = numel(A);
numB = numel(B);
numC = numA + numB;

cntA = ones(numA,1);
cntB = ones(numB,1) + 1;

% sort into a sigle vector

tmp = [A, cntA; B, cntB];

tmp = sortrows(tmp,1)
val = tmp(:,1);
id  = tmp(:,2); % the 1s and 2s

% now move into new vector
C = nan(numC,2); 
c = 1; % counter

for i = 1:numC
   if i == 1
       C(1,id(c)) = val(c);
   else
       if id(i) ~= id(i-1) % stay in same row 
           C(c,id(c)) = val(i);
           c = c + 1; % then advance to next row
       else
           c = c + 1; % advance to next row *first*
           C(c,id(c)) = val(i);
       end
   end
    
end

