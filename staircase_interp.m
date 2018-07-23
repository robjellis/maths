function output = staircase_interp(Yds)

% perform a staircase interpolation of "staircase" data (i.e., inter-step
% intervals that have a low sampling rate ("Yds")
%
% "disc_all" is the series of x- and y-values of the discrete values that
% can then be interpolated
%
% "cont_all" is the interpolated series, which has a default resolution of
% 1000 points (so as to provide better estimates of percentiles)

    ds_thr = .00001; % can adjust this if needed; "0" seems not to work due to small residual decimals
    Yds_sort = sort(Yds);

    % what bin sizes do we *actually* see? 
    Yds_sort_diff = diff(Yds_sort); 
    obs_bins = Yds_sort_diff(Yds_sort_diff > ds_thr); % in theory, all the numbers in obs_bins should be identical


    %% find the index locations of the bin changes
    cc = 1; % counter for this

    for ii = 2:numel(Yds_sort)
        if Yds_sort(ii) - Yds_sort(ii-1) > ds_thr
            bin_ind(cc,1) = ii;
            cc = cc + 1;
        end
    end

    %% now we can create the x and y series to interpolate

    % let's always anchor at the lowest point (and
    % the highest point)
    disc_x(1,1) = 1;
    disc_y(1,1) = Yds_sort(1);

    cc = 2;

    %% horizontal segments
    for ii = 2:numel(bin_ind)
        % *ignore* the first horizontal segment
        disc_x(cc,1) = median(         bin_ind(ii-1) : bin_ind(ii)-1);  % won't necess. be an integer
        disc_y(cc,1) = median(Yds_sort(bin_ind(ii-1) : bin_ind(ii)-1)); % use the median to be more stable 
        cc = cc + 1;
    end

    %% vertical segments; still use the current value of cc and add it to the running list
    for ii = 2:numel(bin_ind)
        % the first vertical is between the first datapoint [disc_y(1)] sand the first change

        disc_x(cc,1) = bin_ind(ii-1) - 0.5;
        disc_y(cc,1) = 0.5 * (disc_y(ii) + disc_y(ii-1)); % use the median values from the *above* step to be stable
        cc = cc + 1;    
    end

    % now we add the highest value
    disc_x(cc) = numel(Yds_sort);
    disc_y(cc) = Yds_sort(end);

    disc_y_tmp = sort(disc_y); % need this for the next step

    % and also the midpoint of the final vertical segment
    disc_x(cc+1) = bin_ind(end)-1;
    disc_y(cc+1) = 0.5 * (disc_y_tmp(end) + disc_y_tmp(end-1));

    %% combine into a matrix
    disc_all = [disc_x disc_y];

    % and SORT to get the original index order
    disc_all = sortrows(disc_all,1);
    
    % check for repeats; RJE should look at the code to see why this
    % happens anyway
    x_disc = disc_all(:,1);
    y_disc = disc_all(:,2);
    
    for i = 2:size(disc_all,1);
       if x_disc(i) == x_disc(i-1)
           % take the average of the y-value, and then delete the previous one
           y_disc(i) = 0.5 * (y_disc(i) + y_disc(i-1));
           x_disc(i-1) = nan;
           y_disc(i-1) = nan;
       end
    end
    
    % then get rid of NaN - careful to include negative values!
    x_disc = x_disc(isnan(x_disc) == 0);
    y_disc = y_disc(isnan(y_disc) == 0);
    
    %% now we do the piecewise cubic hermite interpolating polynomial (PCHIP)
    xi = linspace(1,numel(Yds),1000); % 1000 interpolants to help us resolve percentiles
    
    yiH = pchip(x_disc,y_disc,xi);
    
    %% compare that to regular cubic spline for comparison
    yiS = spline(x_disc,y_disc,xi);
    
    %% outputs
    output.x_disc = x_disc;
    output.y_disc = y_disc;
    output.xi = xi;
    output.yiH = yiH;
    output.yiS = yiS;