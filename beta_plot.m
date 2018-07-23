function beta_plot(a,b)

x = 0:.001:1;
pdf = betapdf(x,a,b);

figure(50)
plot(x,pdf)
xlabel('X')
ylabel('Density')