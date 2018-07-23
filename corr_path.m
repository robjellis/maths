function [result] = corr_path(DATA)

% given a matrix of raw DATA with M rows (subjects) and N columns (measures),
% return the order of columns (values 1:N) which yields the highest
% correlation values
%
% RJE | 2014.05.27

%% first, calculate the correlation; use Spearman by default

corr_mat = corr(DATA,'type','Spearman');
    
    % what is K?
    K = size(DATA,2);
    
    nc_real = NaN(K,K);
    nc_check = NaN(K,K);
    
    
    % the combinations to test
    combs = combnk(1:K,2);
    ncombs = size(combs,1);
    c_col = zeros(ncombs,1);
    r_col = zeros(ncombs,1);
    nc_col = zeros(ncombs,3);
    

    
    for n = 1:ncombs
        c = combs(n,1); % column
        r = combs(n,2); % row
        
        a = M(:,r);
        b = M(:,c);
        
        minab = min(a,b); % rje "check" test: is b identical to a, or active everywhere?
        
        nc_real(r,c) = pearson(a,b); % "pearson" is by rje; doesn't require stats toolbox
        %nc_real(c,r) = nc(r,c); % symmetrical
        
        nc_check(r,c) = pearson(a,minab);
        
        c_col(n) = c;
        r_col(n) = r;
        nc_col(n,1) = c;
        nc_col(n,2) = r;
        nc_col(n,3) = nc_real(r,c);
    end
    

    nc_col = sortrows(nc_col,-3);
    A = nc_col(:,1);
    B = nc_col(:,2);
    AB = [A B];
    NC = nc_col(:,3);
    
    
    % now find the order which makes correlations decrease
    if K == 2 % simple 2-way conjunction
       ord = []; % only two elements, so order doesn't matter
    elseif K > 2
        %nc_temp = nc_col
        kt = AB(1,:); % since we sorted, start with the top
        A(1) = NaN;
        B(1) = NaN;
        AB(1,:) = NaN;
        corr_out(1) = NC(1);

        % see which has the next highest correlation           
        for n = 1:ncombs

            if sum(A(n) == kt)>=1 || sum(B(n) == kt)>=1 % only look at the previous candidates
               indm = n; % find the top-pmost case
               corr_out(2) = NC(indm);
               break % stop the loop
            end

        end

        kt;
        kt2 = [A(indm) B(indm)];
        kt3 = mode([kt kt2]);
        
        %A(indm) = NaN;
        %B(indm) = NaN;
        AB(indm,:) = NaN;
        %NC(indm) = NaN;
        
        % now we have the first three positions
        K_in = 1:K;
        K_out(1) = kt(kt~=kt3);
        K_in = K_in(K_in~=K_out(1));
        K_out(2) = kt3;
        K_in = K_in(K_in~=K_out(2));
        K_out(3) = kt2(kt2~=kt3);
        K_in = K_in(K_in~=K_out(3));
        
        
        %K_out
        %K_in
        % now we know the single path
   
        if K >=3
           kt = K_out(3); % needs to be outside k loop
           for k = 4:K;

                for n = 1:ncombs
                    
                    ab = AB(n,:);
                    cond1 = sum(ab == kt)==1;
                    cond2 = (sum(ab(1) == K_in)==1 || sum(ab(2) == K_in)==1);
                    if cond1 == 1 && cond2 == 1
                       kt = ab(ab~=kt); % critical: redefine kt
                       
                       AB(n,:) = NaN;
                       K_out(k) = kt;
                       K_in = K_in(K_in~=kt);
                       corr_out(k-1) = NC(n);
                       break
                   end
                    
                end
                 
                
           end
           
        end
    end % if K == 2
    nc_col
    result.nc_real = nc_real;
    result.nc_check = nc_check;
    result.nc_col = nc_col;
    result.K_out = K_out;
    result.corr_out = corr_out;
    

