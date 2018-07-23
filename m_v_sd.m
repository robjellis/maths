% mean versus standard deviation in normal distribution

iter = input(' How many iterations?: ');
n    = input(' Number of values per iteration: ');
a    = input(' Mean of random: ');
b    = input(' SD of random: ');


for i = 1:iter
   k = (randn(n,1) + a) * b;
   mk(i) = mean(k);
   sk(i) = std(k);
    
end

figure
plot(mk,sk)
