function [percvol pneg ppos] = percap(vol,upper)
% normalize a distribution based on its 1st and 99th percentiles 
% 
% ndim = 1 means an [A x 1] vector of values
% ndim = 2 means an [A x B] matrix of values
% ndim = 3 means an [A x B x C] matrix of values
%
% will use the spm-based percentile calculator adapted by rje

% will disregard 0-value voxels

% get lower
lower = (100 - upper);

% size of original vol
sized = size(vol);

% figure out dimensionality

if numel(sized) == 2
   if sized(1) == 1 || sized(2) == 1 % either [A x 1] or [1 x A]
      ndim = 1;
   else
      ndim = 2;
   end
elseif numel(sized) == 3
   ndim = 3;
elseif numel(sized) == 4
   ndim = 4;
end

% all **non-zero** values
datavec = vol(:);
datavec = datavec(datavec ~= 0);

% find the percentile bounds

pneg = prctile_spm(datavec,lower); % rje-modified function from SPM
ppos = prctile_spm(datavec,upper);

% first, cap the volume at these values

       % cap the negative
       nvol = pneg * (vol <= pneg);
       
       % cap the positive
       pvol = ppos * (vol >= ppos);
       
       % get the middle
       mid1 = vol > pneg;
       mid2 = vol < ppos;
       midf = mid1 .* mid2;
       
       mvol = vol .* midf;
             
       vol2 = (nvol + mvol + pvol);

% output volume
percvol = zeros(sized);

% perform the operation

if ndim == 1
    
        for a = 1:max(sized)
                 
                val = vol2(a);
                
                if val > 0
                    percvol(a) = val / ppos;
                   
                elseif val < 0
                    percvol(a) = -1 * (val / pneg);
                    
                elseif val == 0
                    % stay at zero
                end    
 
       end % a
    
elseif ndim == 2

    for a = 1:sized(1)
        for b = 1:sized(2)
                
                val = vol2(a,b);
                
                if val > 0
                    percvol(a,b) = val / ppos;
                   
                elseif val < 0
                    percvol(a,b) = -1 * (val / pneg);
                    
                elseif val == 0
                    % stay at zero
                end    
 
        end % b
    end % a
    
    
elseif ndim == 3
    
    for a = 1:sized(1)
        for b = 1:sized(2)
            for c = 1:sized(3)
                
                val = vol2(a,b,c);
                
                if val > 0
                    percvol(a,b,c) = val / ppos;
                   
                elseif val < 0
                    percvol(a,b,c) = -1 * (val / pneg);
                    
                elseif val == 0
                    % stay at zero
                end    
            end % c
        end % b
    end % a

elseif ndim == 4
   
    for a = 1:sized(1)
        for b = 1:sized(2)
            for c = 1:sized(3)
                for d = 1:sized(4)
                
                    val = vol2(a,b,c,d);
                
                    if val > 0
                        percvol(a,b,c,d) = val / ppos;

                    elseif val < 0
                        percvol(a,b,c,d) = -1 * (val / pneg);

                    elseif val == 0
                        % stay at zero
                    end    
                end % d
            end % c
        end % b
    end % a
end

