function yvals = plot_cdf(data,method,looks)


% this is a customized function which will plot the percentage of total
% data elements that are >= x-axis value (method = 'g') or <= x-axis value (method = 'l')


data = data(:); % just to be sure

nlooks = numel(looks);


    
    %looks = linspace(min(data),max(data),nlooks); % values to examine
    
    totdata = numel(data);
    
    nan_to_zero = 1; % do we get rid of these cases, or turn them to zeros?
    
    if nan_to_zero == 1
        data(isnan(data)) = 0;
    else
        data(isnan(data)) = [];
    end
    
    
   
    yvals = zeros(nlooks,1);
    
    for k = 1:nlooks
        if strcmp(method,'g')
        yvals(k) = sum(data >= looks(k)) / totdata; % this is equivalent to cumulative sum
        elseif strcmp(method,'l')
        yvals(k) = sum(data <= looks(k)) / totdata;   
        end
    end
    
    yvals = yvals * 100; % so we get a percentage
   
    
    figure(20)
    
    plot(looks,yvals,'MarkerSize',12,'Marker','.')
    
    % include the x-values
    yvals = [looks' yvals];
    
    num_cases = numel(data)
    