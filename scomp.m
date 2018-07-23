function [sthr] = scomp(type,varargin)

% function [zthr] = pcomp(zdata,zcomp)
%
% one tailed (positive) test

if type == 'n'
zref = varargin{1};
alpha = varargin{2};
sthr = norminv(1 - alpha * (1 - normcdf(zref)));

end



