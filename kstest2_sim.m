function out = kstest2_sim(n1,n2,m1,m2,iter,meth)

sig = nan(iter,1);
%h_mat = zeros(iter,1);


%h_man  = zeros(iter,1);
p_comp = zeros(iter,1);
p_ratio = zeros(iter,1);
p_one   = zeros(iter,1);

stat_all = zeros(iter,1);


tic;
for i = 1:iter
    x = randn(n1,1);
    y = randn(n2,1);
    
    % make sure mean is exact
    x = x - mean(x) + m1;
    y = y - mean(y) + m2;
    
    if strcmp(meth,'orig')
        % don't do anything else
        
    elseif strcmp(meth,'quant')
        minN = min(n1,n2);
        %minN = 50;
        q = 1/minN : 1/minN : 1; % will have minN number of points
        q = q(:);

        %prc = linspace(0,100,minN+1); % the +1 is important
        
        % we need a special way of doing this
        [fx xx] = ecdf(x);
        [fy xy] = ecdf(y);
        
        qx = nan(minN,1);
        qy = nan(minN,1);
        
        % note: there may be subtle decimal differences that we want to
        % account for; just round to 3rd decimal place
        mult = 10000;
        q = round(q*mult)/mult;
        fx = round(fx*mult)/mult;
        fy = round(fy*mult)/mult;
        
        for j = 1:minN
          qx(j) = max(xx(fx <= q(j))); 
          
          qy(j) = max(xy(fy <= q(j)));  
        end
        
        %sum(q(:) - fx(2:end))
        
        %zz = [sort(x) qx]
        
        % replace the vars
        x = qx;
        y = qy;
        
    end

    % get the means as a summary statistic
    mx = mean(x);
    my = mean(y);
    
    [h p stat] = kstest2(x,y,.05,'unequal'); % two-tailed p;

    stat_all(i) = stat;
    p_comp(i,1) = p;

    % -------------------------------------------
    % copied from Matlab
%    n1     =  length(x1);
%    n2     =  length(x2);
    n      =  n1 * n2 /(n1 + n2);
    lambda =  max((sqrt(n) + 0.12 + 0.11/sqrt(n)) * stat , 0);
    
    p_one(i) =  exp(-2 * lambda * lambda);
    
    % asymptotic Q function
    j       =  (1:101)';
    pValue  =  2 * sum((-1).^(j-1).*exp(-2*lambda*lambda*j.^2));
    pValue  =  min(max(pValue, 0), 1);  
    p_two_tail = pValue;
    

    p_ratio(i) =  p_two_tail / p_one(i);

    p_comp(i,2) = pValue;
    
%     % RJE manual checking
%     D_crit = 1.36 * sqrt((n1+n2)/(n1*n2));
%     
%     if stat > D_crit
%         h_man(i) = 1;
%     end
    
% revised
    if p < .05 && my > mx
        sig(i) = 1;   % 0.5;
    elseif p < .05 && my < mx
        sig(i) = -1;  %-0.5;
    else
        sig(i) = 0;
    end
    
end % iter loop
tocc = toc;


if iter == 1
    figure(100)
    ecdf(x)
    hold on
    ecdf(y)
    hold off
    ylabel('Cumulative prob.')
    xlim([-4 4]) 
    
    size(x)
    size(y)
end


sig_prc  = 100 *sum(sig) / iter;
%sig_prc = 50 + (100 * sum(sig) / iter);

% h_mat_sig_prc = 100 * sum(h_mat) / iter;
% h_man_sig_prc = 100 * sum(h_man) / iter;

out.sig_prc = sig_prc;
%out.h_mat_sig_prc = h_mat_sig_prc;
%out.h_man_sig_prc = h_man_sig_prc;
out.duration = tocc;
out.p_comp = p_comp;

figure(10)
plot(p_one,p_ratio,'.')

figure(11)
ecdf(stat_all)
xlim([0 1])