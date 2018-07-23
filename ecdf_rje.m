function [Fout,x] = ecdf_rje(y,method)

% Empirical (Kaplan-Meier) cumulative distribution function.
%   [F,X] = ECDF(Y) calculates the Kaplan-Meier estimate of the
%   cumulative distribution function (cdf), also known as the empirical
%   cdf.  Y is a vector of data values.  F is a vector of values of the
%   empirical cdf evaluated at X.
%
%   ECDF(...) without output arguments produces a plot of the empirical
%   cdf. 
%
% Adapted from Matlab's ecdf.m by RJE | 2017.02.11

if nargin < 2
    method = 'cdf';
end

% make sure a vector
x = y(:);
freq = ones(size(x));

% Remove missing observations indicated by NaN's.
t = ~isnan(x);

x = x(t);
n = length(x);
if n == 0
    error('stats:ecdf:NotEnoughData',...
          'Input sample has no valid data (all missing values).');
end

freq = freq(t);

% Sort observation data in ascending order.
[x,t] = sort(x);
%cens = cens(t);
freq = freq(t);

% Compute cumulative sum of frequencies
totcumfreq = cumsum(freq);
obscumfreq = totcumfreq;

t = (diff(x) == 0);
if any(t)
    x(t) = [];
    totcumfreq(t) = [];
    obscumfreq(t) = [];
end
totalcount = totcumfreq(end);

% Get number of deaths and number at risk at each unique X
D = [obscumfreq(1); diff(obscumfreq)];
N = totalcount - [0; totcumfreq(1:end-1)];

% No change in function except at a death, so remove other points
t = (D>0);
x = x(t);
D = D(t);
N = N(t);

if strcmp(method,'cdf') || strcmp(method,'survivor')
    % Use the product-limit (Kaplan-Meier) estimate of the survivor
    % function, transform to the CDF.
    S = cumprod(1 - D./N);
    if strcmp(method,'cdf')
        Func = 1 - S;
        F0 = 0;       % starting value of this function (at x=-Inf)
        funcdisplayname = 'F(x)';
    elseif strcmp(method,'survivor')
        Func = S;
        F0 = 1;
        funcdisplayname = 'S(x)';
    end
elseif strcmp(method,'cumhazard')
    % Use the Nelson-Aalen estimate of the cumulative hazard function.
    Func = cumsum(D./N);
    F0 = 0;
    funcdisplayname = 'H(x)';
end

% Include a starting value; required for accurate staircase plot
x = [min(y); x];
F = [F0; Func];

%% Plot if no return values are requested
if nargout == 0
    h = stairs(x,F);
    xlabel('x');
    ylabel(funcdisplayname);
    %Set custom data cursor text for data line.
    hB = hggetbehavior(h(1),'datacursor');
    set(hB,'UpdateFcn',@ecdfDatatipCallback);
    setappdata(h(1),'funcdisplayname',funcdisplayname);
    setappdata(h(1),'D',[0; D]);
    
else
    Fout = F;
end

