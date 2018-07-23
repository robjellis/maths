% to see how many data points makes a stable estimate of SD


iter = input(' How many interations?');
points = 100;

stdev = zeros(iter,points-1);

for i = 1:iter
    
   a=rand(points,1);
   
   for j = 1:points-1
       stdev(i,j) = std(a(1:j+1));
   end
    
    
end

figure
hold on
for i=1:iter
plot(1:j,stdev(i,1:j))
fin_std = std(stdev);
end
hold off

figure
plot(fin_std)