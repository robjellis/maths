%
% make histogram outlines for N t distributions with v d.f.

df = input(' Enter the d.f. for the T-distribution: ');
N = input(' Enter the size of each T-distribution: ');
j = input(' Enter the number of simulations to run: ');

x = -6:.25:6;  % x values for bins

figure
hold on
for i = 1:j
    
t = trnd(df,N,1);
y = histc(t,x);
plot(x,y);
hold on

end


fprintf('\n Finished.\n\n');