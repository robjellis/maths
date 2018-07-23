% bootstrapping to compare observed r-values between a fixed x-variable
% (IV) and random permutations of a y-variable


numx = numel(x);
numy = numel(y);

if numx ~= numy
    return
end

iter = input(' How many iterations to perform?: ');


rvals = zeros(1,iter);

for i = 1:iter
   newy = zeros(numy,1);
   yrand = randperm(numy);
   for j = 1:numy
       newy(j) = y(yrand(j));
   end
   rvals(i) = corr(x,newy);
   
end