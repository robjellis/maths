function output = ranksum_sim(nfull,nsamp,dist,iter1,iter2)

% let's look at how distribution properties affect the Z-test of ranksum

full1_mn = nan(iter1,1);
full1_sd = nan(iter1,1);
Z_md     = nan(iter1,1);
T_md     = nan(iter1,1);

progressbar(0,0)

for i = 1:iter1
    
    if strcmp(dist,'u')
        c = rand(1,1); % random from uniform dist
        d = rand(1,1);

        a = min(c,d);
        b = max(c,d);

        full1 = a + (b-a).*rand(nfull,1);
        full1 = full1*100;

    elseif strcmp(dist,'ureal')
        full1 = rand(nfull,1) + rand;
        full1 = full1*100;
        
    elseif strcmp(dist,'n')
        full1 = randn(nfull,1) + rand;

    end
    
    full1_mn(i) = mean(full1);
    full1_sd(i) =  std(full1);
    
    % subsampling
    Z_temp = nan(iter2,1);
    
    T_temp = nan(iter2,1);
    
    for j = 1:iter2
        
        ind1  = randperm(nfull); % all values
        ind1 = ind1(1:nsamp); % just the first ones
        samp1 = full1(ind1); 


        if strcmp(dist,'u') || strcmp(dist,'ureal')
            % sample 2 is from the regular distribution
            samp2 = rand(nsamp,1); % just make this large
            samp2 = samp2*100;

        elseif strcmp(dist,'n')
            samp2 = randn(nsamp,1);
            
        end
        % do the MWW test
        [p h stats] = ranksum(samp1,samp2);
        Z_temp(j) = stats.zval;
        
        % two sample T-test
        [h p ci stats] = ttest2(samp1,samp2);
        T_temp(j) = stats.tstat;

        progressbar([],j/iter2)
    end
    
    Z_md(i) = median(Z_temp);
    T_md(i) = median(T_temp);
    
    progressbar(i/iter1,[])
end
progressbar(1,1)

figure(200)
plot(full1_mn,Z_md,'.')
xlabel('Samp1 mean')
ylabel('Samp1 vs Samp2 MWW Z-value')

if strcmp(dist,'u')
   xlim([0 100]) 
   ylim([-5 5])
end

figure(201)
plot(full1_mn,T_md,'.')
xlabel('Samp1 mean')
ylabel('Samp1 vs Samp2 T-test')

if strcmp(dist,'u')
   xlim([0 100]) 
end

% fit a line
coef = polyfit(full1_mn,Z_md,1);
rho1 = corr(full1_mn,Z_md,'type','Spearman');
rho2 = corr(full1_mn,T_md,'type','Spearman');

output.coef_MWW = coef;
output.rho_MWW = rho1;
output.rho_T   = rho2;
