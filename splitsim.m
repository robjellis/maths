function [stats] = splitsim(cat,samp,iter,Amult,Bmult,figs)
% 
% function [res] = splitsim(cat,samp,iter,Amult,Bmult,figs)
%
% simulation to show how power is lost when doing an n-way split of a continuous variable
%
% cat = input(' Number of categories: ');
% samp = input(' Number of cases per category: ');
% iter = input(' Number of iterations: ');
% mult = input(' Multiplier for y values: ');
%
% RJ Ellis, Nov 2011


r1 = zeros(iter,1);
r2 = zeros(iter,1);
z1 = zeros(iter,1);
z2 = zeros(iter,1);

num = cat * samp;
x2 = zeros(num,1);

k = 1;
for j = 1:cat
    
x2(k:k+samp-1) = j;

k = k + samp;
   
end

for i = 1:iter
    
x1 = randn(num,1);  % original data
x1 = sort(x1);      % sort has to happen to correspond with sorted categorical vector

y = x1*Amult + randn(num,1)*Bmult;





r1(i) = corr(x1,y);
r2(i) = corr(x2,y);

%bins = -1:.025:1;
%p1 = histc(r1,bins);
%p2 = histc(r2,bins);

% now convert r-values to z-values so we get a truer sense of things
z1t = corr(x1,y);
z2t = corr(x2,y);
z1(i) = 1/2 * log((1+z1t)/(1-z1t));
z2(i) = 1/2 * log((1+z2t)/(1-z2t));
end


bins = -1:.025:3;
p1 = histc(z1,bins);
p2 = histc(z2,bins);

figure(70)
subplot(3,1,1)
if figs == 1
    plot(x1,y,'.');
    axis([-3 3 -3 3])
elseif figs == 0
    aa=0;
    bb=0;
    plot(aa,bb);
end
subplot(3,1,2)
plot(bins,p1,'g') % the proper correlation
hold on
plot(bins,p2,'r') % the categorical split
hold off

subplot(3,1,3)
plot(z1,z1-z2,'.')

% statistics

stats.mean_cont_r = mean(r1);
stats.mean_split_r = mean(r2);
stats.mean_cont_z = mean(z1);
stats.mean_split_z = mean(z2);
