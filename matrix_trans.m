%for spatially aligning matrices of different sizes

c=zeros(9,9,9);  % a large empty matrix with origin at 5,5,5
centc = [5 5 5]; % origin of the empty matrix

a = rand(3,5,7);   % center is at 2,3,4;
centa = [2 3 4];
[d1 d2 d3] = size(a);
tx = centc(1)-centa(1); 
ty = centc(2)-centa(2);
tz = centc(3)-centa(3);
  

for xi=1:d1
    for yi=1:d2
        for zi=1:d3
            c(tx + xi, ty + yi, tz + zi) = a(xi,yi,zi);
        end
    end
end

