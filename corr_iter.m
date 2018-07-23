% correlation random distributions


vals = input(' How many values in each vector?: ');
iter = input(' How many iterations?: ');

r = zeros(iter,1);

for i = 1:iter
    
    x = randn(iter,1);
    y = randn(iter,1);
    
    r(i) = corr(x,y); % the r value
end

hist(r,50);
