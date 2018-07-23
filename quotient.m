function output = quotient(dtype,M,N,iter)

% to try out the new method by rje
%
% quotient(dtype,M,N,iter)

prc = 50:10:90;

all_loc_prc = nan(iter,numel(prc));
all_suc_prc = nan(iter,numel(prc));

for j = 1:iter
    
    if strcmp(dtype,'n')
        x = randn(N,1)*20+M*1000;
        x = x/1000;

    elseif strcmp(dtype,'e')
        x = exprnd(1,N,1);
        x = x + M;

    end
    x = x(:);

    % outlier

    %x(100) = max(x) * 3;

    % quotient versus midhinge
    loc = 0.5 * (prctile_nist(x,75) + prctile_nist(x,25));


    % 3. ratio deviation from loc

    for i = 1:numel(x)
        if x(i) >= loc
            quot_loc(i) = x(i) / loc;
        else
            quot_loc(i) = loc / x(i);
        end
    end

    quot_loc(1) = NaN;
    quot_loc = quot_loc(:);


    % new SUCCESSIVE quotient method

    for i = 2:numel(x)
        if x(i) >= x(i-1)
            quot_suc(i) = x(i) / x(i-1);
        else
            quot_suc(i) = x(i-1) / x(i);
        end
    end

    quot_suc(1) = NaN;
    quot_suc = quot_suc(:);


%% alternative: using log

% take the log
log_x = log(x);

% take the diff
diff_log_x = diff(log_x);

% get the median
med = median(diff_log_x); 

    log_suc = exp(abs(diff(log(x)))); % identical to above method!
    

    
    
%% plots
    if iter < 10
        figure(500)
        subplot(5,1,1)
        plot(x);    
        ylabel('x')
        
        subplot(5,1,2)
        plot(sort(diff(x)))
        ylabel('sort(diff((x))')

        subplot(5,1,3)
        plot(quot_suc);
        ylabel('succ. quot.')

     
        
        figure(501)
        % log version
        
        subplot(5,1,1)
        plot(x);
        ylabel('x')
        
        subplot(5,1,2)
        plot(sort(diff_log_x))
        ylabel('sort(log(x))')
        
        subplot(5,1,3)
        plot(log_suc)
        ylabel('log_suc')
        
        figure(502)
        plot(quot_suc(2:end),log_suc,'.')

        % correlation
        figure(510)
        plot(quot_loc(2:end),quot_suc(2:end),'.')
        corr_meas = corr(quot_loc(2:end),quot_suc(2:end));
    end

    % now percentiles

    quot_loc_prc = prctile_nist(quot_loc(2:end),prc);
    quot_suc_prc = prctile_nist(quot_suc(2:end),prc);

    all_loc_prc(j,:) = quot_loc_prc';
    all_suc_prc(j,:) = quot_suc_prc';
    
    % get the equiv of std

    % x_dev = median(x) * (xprc - 1)
    % 
    % x_std = std(x)

    % now multiply by x

    %x_final = median(x) * xprc;
   

end

% correlations
figure(510)

for p = 1:numel(prc)
   subplot(1,numel(prc),p)
   plot(all_loc_prc(:,p),all_suc_prc(:,p),'.')
   
   all_corr(1,p) = corr(all_loc_prc(:,p),all_suc_prc(:,p));
end

%% outputs


if iter < 10
    output.quot_loc_prc = quot_loc_prc';
    output.quot_suc_prc = quot_suc_prc';
    output.corr = corr_meas;
else
    output.corr = all_corr;
end
