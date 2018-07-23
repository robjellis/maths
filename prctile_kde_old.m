function output = prctile_kde_old(data,kde_meth,prc_val,plot_it)

%
% "data" assumes a single vector of numbers
%
% "prc_val" is 0 to 100 scale
%

if size(data,2) > 1
    % we only want a vector
    fprintf('\n Warning: this function only accepts an [N x 1] vector of data.\n\n')
    return
end

if nargin < 3
    prc_val = 1:99;
end

if nargin < 4
    plot_it = 0;
end

data = data(:);

%% KDE estimation
% http://www.mathworks.com/matlabcentral/fileexchange/14034-kernel-density-estimator

nmesh = 2^12; % default = 2^12

xinds = 1:nmesh;

% we need a COMMON scale for the data
range = max(data) - min(data);
MIN = min(data) - range/4; 
MAX = max(data) + range/4;

if strcmp(kde_meth,'m')
    % matlab
    X = linspace(MIN,MAX,nmesh);

    [D]   = ksdensity(data,X,'function','pdf');
    [CDF] = ksdensity(data,X,'function','cdf');

elseif strcmp(kde_meth,'b')
    % botev
    [bw,D,X,CDF] = kde(data,nmesh,MIN,MAX); % control
end

% make sure we are bounded correctly
D(D < 0) = 0;
CDF(CDF < 0) = 0;
CDF(CDF > 1) = 1;

%% percentiles

% standard one
output.nist_percentiles = prctile_nist(data,prc_val);

p_val = prc_val/100; % 0 to 100 -> 0 to 1
nump = numel(p_val);

for p = 1:nump
   prc(p) = min(X(CDF >= p_val(p))); % the x-value
end

output.kde_percentiles = prc;


%% plots

if plot_it == 1
   figure(500)
   subplot(1,4,1)
   nbins = ceil(numel(data)/10);
   if nbins > 30
      nbins = 30;
   end
   hist(data,nbins)
   xlabel('Data')
   ylabel('Count')
  
   y = hist(data,nbins);
   axis([MIN MAX 0 max(y)*1.1])

   % PDF
   subplot(1,4,2)
   plot(X,D,'b')
   xlabel('Data')
   ylabel('Estimated probability')
   title('PDF')
   ymax = max(max(D));
   axis([MIN MAX 0 ymax*1.1])

   % CDF
   subplot(1,4,3)
   plot(X,CDF,'b')
   xlabel('Data')
   ylabel('Pr(D<=X)')
   title('CDF')
   axis([MIN MAX 0 1]) 
   
   % ECDF
   subplot(1,4,4)
   ecdf(data)

end
