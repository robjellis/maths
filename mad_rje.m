function output = mad_rje(x)

% Mean or Median absolute deviation, with adjustment to match scale of
% standard deviation ()
%
% * standard formula, but included here to avoid using the Matlab Stats toolbox
% * for simplicity, we only perform calculation on a single vector ([N x 1])

x = x(:);

% get rid of NaN if present
x(isnan(x) == 1) = [];

% check size
if isempty(x) == 0
    % OK

    % Compute the median of the absolute deviations from the median.
    mdad = median(abs(x - median(x))) * 1.4826; % the extra term is to get the scale appropriate to match STDEV 
                                                % http://en.wikipedia.org/wiki/Median_absolute_deviation#Relation_to_standard_deviation
                                                % 1.4826 = 1/norminv(0.75)

    % Compute the mean of the absolute deviations from the mean.
    mnad = mean(abs(x - mean(x))) * sqrt(pi/2); % the extra term is to get the scale appropriate to match STDEV 
                                                % http://en.wikipedia.org/wiki/Absolute_deviation

    output.std = std(x);
    output.mean = mean(x);
    output.mnad = mnad;
    output.cvmnad = 100 * mnad / mean(x);
    output.median = median(x);
    output.mdad = mdad;
    output.cvmdad = 100 * mdad / median(x);
    
else
    %fprintf(' Warning: this function only operates on a single [N x 1]) vector. \n');
    output.std = nan;
    output.mean = nan;
    output.mnad = nan;
    output.cvmnad = nan;
    output.median = nan;
    output.mdad = nan;
    output.cvmdad = nan;  
end