function [inv] = percit(data,ndim)
% normalize a distribution based on its 1st and 99th percentiles 
% 
% ndim = 1 means an A x 1 vector of values
% ndim = 2 means an A x B matrix of values
% ndim = 3 means an A x B x C matrix of values
%
% will disregard 0-value voxels

% size of original data
sized = size(data);

% get the data into a column
datavec = data(:);

% find the 1%ile and 99%ile

p1 = prctile(datavec,1);
p99 = prctile(datavec,99);

% just the pos values
datapos = datavec(datavec > 0);

% just the negative values
dataneg = datavec(datavec < 0);

% how many positive and negative voxels?
numpos = numel(datapos);
numneg = numel(dataneg);

% output volume
inv = zeros(sized);

% get the inverse percentiles

if ndim == 1
    
    
elseif ndim == 2

    for a = 1:sized(1)
        for b = 1:sized(2)
        end
    end
    
    
elseif ndim == 3
    
    for a = 1:sized(1)
        for b = 1:sized(2)
            for c = 1:sized(3)
                
                val = data(a,b,c);
                
                if val > 0
                    inv(a,b,c) = sum(datapos < val) / numpos;
                   
                elseif val < 0
                    inv(a,b,c) = -1 * (sum(dataneg > val) / numneg);
                    
                elseif val == 0
                    % stay at zero
                end    
            end % c
        end % b
    end % a

end

