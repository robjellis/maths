function Xds = downsample_rje(X,type,Hz,plot_it)

% downsample the vector X, which can either be 
%     1. inter-event intervals (type = 'i'), in *seconds* or 
%     2. timestamps (type = 't'), in *seconds*
%
% 'Hz' is the value of the desired downsampled series (e.g., 'Hz'= 50)
%
% 'Xds' will have the same format (either 'i' or 't' as the input format)
%
% rje | 2013.10.02

X = X(:);

% convert X to timestamps if needed

if strcmp(type,'i')
    Xts = cumsum(X);
elseif strcmp(type,'t')
    % do nothing since we already have timestamps
    Xts = X;
end

% convert Hz into sampling rate
sr = (1 / Hz); % samples every 'sr' *seconds*

ds = ceil(Xts / sr) * sr; % the new time stamp is rounded *up* to the next sampling point

% and now we convert back to the input format

Xds = nan(size(X)); % initialize this variable

if strcmp(type,'i')
   Xds(1) = ds(1);
   Xds(2:numel(ds)) = diff(ds); % this is now the resampled ITI series, which is compared to original X series
   
elseif strcmp(type,'t')
    % do nothing since we keep format in timestamps
end

if plot_it == 1
   figure(100)
   plot(X,'b','LineWidth',2)
   hold on
   plot(Xds,'r')
   hold off
end
