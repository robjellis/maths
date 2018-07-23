function [pval] = prctile_inv(data,xtar,tail)

% function [pval] = prctile_inv(data,xtar,tail)
%
% calculate the "inverse percentile": given a distribution of values,
% return the corresponding percentile (i.e., proportion of values < x) for
% a target value (xtar)
%
% of course, this only becomes useful if DATA is large (e.g., > 10000 values)
%
% tail is 'l' 'r' or 'b'
%
% rje, july 2012

data = data(:); % just in case

% how many targets?
xtar = xtar(:);
ntar = numel(xtar);

pval = zeros(ntar,1);

% let's get the percentile values

prcs = (0:.01:100)'; % will give p-values in steps of .0001

xobs = prctile(data,prcs);

% now just find the value of xobs closest to xtar

for i = 1:ntar
    idx = max(find(xobs<=xtar(i)));

    % now get the prctile value for this
    prc = prcs(idx);


    if strcmp(tail,'l')
        pval(i) = prc/100;

    elseif strcmp(tail,'r')
        pval(i) = 1 - prc/100;

    elseif strcmp(tail,'b')
        if prc < 50
            pval(i) = prc/100 * 2;
        elseif prc >= 50
            pval(i) = (1 - prc/100) * 2;    
        end
    end

end


