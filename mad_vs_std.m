function output = mad_vs_std(iter)

% visualizing the relationship as a function of sample size

num = [8 16 32 64 128 250 512 1024 2048 4096 8192 16384];

%% loop

% variables

res = nan(numel(num),4);

alls = nan(numel(num),4);
allm = nan(numel(num),4);

for n = 1:numel(num)    
    x = randn(num(n),iter); 
    s = std(x); 
    m = mad(x,1); 
    
    rat = s ./ m;
    rat = rat(:);
    
    prc = [10 50 90];
    perc = prctile_nist(rat,prc);
    
    percs = prctile_nist(s(:),prc);
    percm = prctile_nist(m(:),prc);   
    
    res(n,1) = num(n);
    res(n,2) = perc(1);
    res(n,3) = perc(2);
    res(n,4) = perc(3);
    
    % save the STD values
    alls(n,1) = num(n);
    alls(n,2) = percs(1);
    alls(n,3) = percs(2);
    alls(n,4) = percs(3);   
    
    % save the MAD values
    allm(n,1) = num(n);
    allm(n,2) = percm(1);
    allm(n,3) = percm(2);
    allm(n,4) = percm(3);     
    
end

%% plot

%% output
output.ratio = res;
output.std = alls;
output.mad = allm;