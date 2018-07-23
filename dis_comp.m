% 3-D Euclidian distance (for brain coordinates)
%
% Batch use: will prompt for first set and second set
%
% rje 2012

function [d] = dis_comp

intyp = input('\n Enter coordinates via [1] copy-paste or [2] text files: ');
if intyp == 1
    a = input('\n Enter 1st set of XYZ [N x 3]: ');
    b = input(' Enter 2nd set of XYZ [N x 3]: ');
elseif intyp == 2
   f = spm_select(2,'mat','Select the Batch 1 and Batch 2 coordinate files:',[],pwd,'.*');
   fa = f(1,:); 
   fa = strtrim(fa);
   a = load(fa);
   fb = f(2,:); 
   fb = strtrim(fb);
   b = load(fb);
    
end
numa = size(a,1);
numb = size(b,1);

if numa == numb
    % OK
else
    fprintf(' \n Warning: Number of coordinates in 1st and 2nd sets does not match.\n');
    return
end

d = zeros(numa,1);

for i = 1:numa
    
    x1 = a(i,1); y1 = a(i,2); z1 = a(i,3);
    x2 = b(i,1); y2 = b(i,2); z2 = b(i,3);

    d(i) = sqrt((x2 - x1)^2 + (y2 - y1)^2 + (z2 -z1)^2);
    
end