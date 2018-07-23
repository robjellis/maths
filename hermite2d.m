function [output] = hermite2d(x,data,xint)

% for a R x C data matrix, output will return the pchip results by looping
% through each row
%
% x should be a 1 x C vector
% int can be any size, with the same range as x

% check to make sure size is OK

sized = size(data);
sizex = size(x);
sizei = size(xint);

if sized(2) == sizex(2)
    % OK
else
    fprintf('\n Warning: the data and the interpolant do not have the same number of columns.\n');
    return
end

output = zeros(sized(1),sizei(2));

for i = 1:sized(1)
    datavec = data(i,:);
    output(i,:) = pchip(x,datavec,xint);
    
end


