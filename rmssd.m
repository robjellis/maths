function output = rmssd(x)
% 
% x = a vector of IBIs (*not* raw timestamps)
%
% ** note: function td_calcs is a more recent version of this!
%
% measures calculated:
% - CVSD
% - CVRMS
% - PAD: percentile absolute deviation from median IBI
% - PSD: percentile successive difference from the median successive difference
%
% rje | version = 2013.03.01

if isnan(x) 
     RMS = NaN;
   CVRMS = NaN;
   SD = NaN;
   CV = NaN;
else  

       RMS = sqrt(mean(diff(x).^2));  % find all pairwise differences, square them, take the mean, and take the square root
     CVRMS = (RMS / mean(x)) * 100;  % rms relative to the mean (similar to coefficient of variation)
    
        SD = std(x);
        CV = SD / mean(x) * 100; % note: if mean(x) approaches zero, then CV will be highly inflated
        

end

output.rms = RMS;
output.cvrms = CVRMS;
output.sd = SD;
output.cvsd = CV;