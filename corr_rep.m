function [stats] = corr_rep(vec,N,rootit,ceilit)

% theoretical z values for a "repeated measures correlation"
% function [stats] = corr_rep(vec,N)

% input vector
vec = vec(:);

num = numel(vec);

data = randn(num,N);

% do the correlation at all N
r = corr(vec,data);

% fisher transform

r2z = 0.5*log((1+r)./(1-r));

% optional: take square root abs value of z-vals > 1
% rje: this doesn't produce a nice distribution shape!
if rootit == 1
    snz = sign(r2z);
    abz = abs(r2z);
    gt1 = abz .* (abz > 1);
    lt1 = abz .* (abz <= 1);
    
    sqz = sqrt(gt1);
    
    r2z = (sqz + lt1) .* snz;
    
else
    r2z = r2z;
    
end

if ceilit == 1
   snz = sign(r2z);
   abz = abs(r2z);
   gt5 = 5 * (abz > 5); % caps it at 5
   lt5 = abz .* (abz <= 5);
   
   z = (lt5 + gt5) .* snz;
   
else
   z = r2z;
end

% plots

figure(80)
hist(r,100)

figure(81)
hist(z,100)


% stats

stats.min = min(z);
stats.mean = mean(z);
stats.max = max(z);
stats.std = std(z);
stats.skew = skewness(z);
stats.prctile_95 = prctile(z,95);
stats.prctile_99 = prctile(z,99);



