function output = prctile_kde(data,kde_meth,bw,prc_val,range,plot_it)

%
% "data" assumes a single vector of numbers
%
% "prc_val" is the set of target percentile values (e.g., 1 to 99)
%

if size(data,2) > 1
    % we only want a vector
    fprintf('\n Warning: this function only accepts an [N x 1] vector of data.\n\n')
    return
end


if nargin < 3
    prc_val = 1:99;
end

if isempty(range)
    % make the minimum zero if all values are 0 or positive
    if sum(data < 0) == 0
        MIN = 0;
        ext = (max(data) - 0)/2;
        MAX = max(data)+ext;
    else
        ext = (max(data) - min(data)) / 2; % div by 10 is standard per kde.m; RJE made it wider
        MIN = min(data) - ext;
        MAX = max(data) + ext;
    end

else
    MIN = range(1);
    MAX = range(2);
end

nump = numel(prc_val);

if nargin < 6
    plot_it = 1;
end

data = data(:);

%% KDE estimation
% http://www.mathworks.com/matlabcentral/fileexchange/14034-kernel-density-estimator



if strcmp(kde_meth,'m')
    % matlab
    nmesh = round(numel(data)/10); % 100 is default but this is too coarse
    X = linspace(MIN,MAX,nmesh);
    
%     X = prc_val/100; % goes from 0 to 1 
%     X = [X(:); 1];
%     X = unique(X);
    
    if isempty(bw)
        [D zzz bw]   = ksdensity(data,X,'support',[MIN MAX],'npoints',nmesh,'function','pdf');
        [CDF zzz bw] = ksdensity(data,X,'support',[MIN MAX],'npoints',nmesh,'function','cdf');          
    else
        [D zzz bw]   = ksdensity(data,X,'support',[MIN MAX],'width',bw,'npoints',nmesh,'function','pdf');
        [CDF zzz bw] = ksdensity(data,X,'support',[MIN MAX],'width',bw,'npoints',nmesh,'function','cdf');    
    end


elseif strcmp(kde_meth,'b')
    % botev
    nmesh = 2^12; % default per Botev
    [bw,D,X,CDF] = kde(data,nmesh,MIN,MAX); % control
end

% make sure we are bounded correctly

% find the first point where CDF hits 1.0, and make sure all subsequent
% values are 1.0

if min(CDF) < 0
    this_ind = find(CDF < 0,1,'last'); % find the furthest value
    CDF(1:this_ind) = 0;    
else
    % don't do anything
end

if max(CDF) > 1
    this_ind = find(CDF > 1,1,'first'); % find the earliest one
    CDF(this_ind:end) = 1;
else
    % don't have to do anything
end

%% percentiles

% standard one

prc_simple = nan(nump,1);
prc_kde    = nan(nump,1);

for p = 1:nump
   this_p = prc_val(p); 
   prc_simple(p) = prctile(data,this_p);

   prc_kde(p) = min(X(CDF*100 >= this_p)); % the x-value
end

%% plots

if plot_it == 1
   figure(500)
   subplot(1,2,1)
   nbins = ceil(numel(data)/10);
   if nbins > 50
      nbins = 50;
   end
   hist(data,nbins)
   xlabel('Data')
   ylabel('Count')
   title('Histogram')
  
   y = hist(data,nbins);
   axis([MIN MAX 0 max(y)*1.1])

   % PDF
   subplot(1,2,2)
   plot(X,D,'b')
   xlabel('Data')
   ylabel('Probability')
   title('kde PDF')
   ymax = max(max(D));
   axis([MIN MAX 0 ymax*1.1])

   figure(502)
   % CDF
   subplot(1,2,1)
   plot(X,CDF,'b')
   xlabel('Data')
   ylabel('Cumulative prob.')
   title('kde CDF')
   axis([MIN MAX 0 1]) 
   
   % ECDF
   subplot(1,2,2)
   ecdf(data)
   title('Data ECDF')
   axis([MIN MAX 0 1]) 
   
   figure(502)
   plot(prc_simple,prc_kde)
   mm = min(min(prc_simple),min(prc_kde));
   nn = max(max(prc_simple),max(prc_kde));
   axis([mm nn mm nn])
   xlabel('Simple percentile')
   ylabel('KDE percentile')
end

%% outputs
output.bandwidth = bw;
output.prc_target = prc_val(:);
output.prc_simple = prc_simple;
output.prc_kde    = prc_kde;
