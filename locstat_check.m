function locstat_check(iter)

% 

N = [20 40 60 80 100 250 500 1000 5000];

numn = numel(N);

pdiff_all = nan(iter,numn);

for i = 1:iter
    
    for j = 1:numn
        x = randn(N(j),1);

        out = locstat(x);

        pdiff_all(i,j) = out.pdiff;
    end
    
end

figure(10)
plot(N,prctile_nist(pdiff_all,95),'MarkerSize',10,'Marker','.')
xlabel('Number of data points')
ylabel('95th percentile of percent difference')

% boxplot_rje(pdiff_all,[1 99],10)
    
    