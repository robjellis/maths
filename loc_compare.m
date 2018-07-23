function output = loc_compare(X,plot_it)

% compare the various methods of getting the location

%% standard time domain

loc_mean = mean(X);
loc_median = median(X);
loc_mode = mode(X);
loc_min = min(X);

loc_trimean = (prctile_nist(X,25) + 2 * prctile_nist(X,50) + prctile_nist(X,75))/4;

%% KDE (two methods)

npoints = 2^10; % used for both methods

% Matlab KDE
tic;
        range = max(X) - min(X);
        MIN = min(X) - range/4; 
        MAX = max(X) + range/4;

        [f,xi] = ksdensity(X,linspace(MIN,MAX,npoints)); 

        % simple rule: if we don't have enough unique, then just use mode
        prc_unique = 100 * numel(unique(X)) / numel(X);
        
        if prc_unique < 50
            % just use mode to be safe
            kde_matlab = mode (X);
        else
            kde_matlab = min(xi(f == max(f))); % take min so we get the faster rate
        end
    
kde_matlab_time = toc;

% Botev
tic;
[bandwidth,density,xmesh]= kde(X,npoints); % use same number of points as Matlab
kde_botev = min(xmesh(density == max(density)));
kde_botev_time = toc;

%% RJE method
plot_brute = 0;
%brute_out = cent_tend_brute(X,1,0,plot_brute);
%loc_brute = brute_out.brute_location;

%% plots

if plot_it == 1
    figure(600)
    nplot = 5;

    
    subplot(1,nplot,1)
    plot(X)
    title('Original')

    subplot(1,nplot,2)
    plot(sort(X))
    title('Sorted')

    %subplot(1,nplot,3)
    %ecdf(X)

    subplot(1,nplot,3)
    hist(X,round(numel(X)/10))
    title('Histogram')

    subplot(1,nplot,4)
    plot(xi,f);
    title('Matlab KDE')

    subplot(1,nplot,5)
    plot(xmesh,density)
    title('Botev KDE')
end

%% outputs

output.mean = loc_mean;
output.trimean = loc_trimean; 
output.median = loc_median;
output.mode = loc_mode;
%output.loc_brute = loc_brute;
output.kde_matlab = kde_matlab;
output.kde_matlab_time = kde_matlab_time;
output.kde_botev = kde_botev;
output.kde_botev_time = kde_botev_time;
output.min = loc_min;
output.percent_unique = 100*numel(unique(X))/numel(X);

