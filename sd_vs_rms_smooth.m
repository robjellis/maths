function [output] = sd_vs_rms_smooth(iter,N,M,CV,alpha,Hz)

% compare sd (CV) with rms (CV)
%
% [output] = sd_vs_rms(iter,N,M,CV,alpha,Hz)

cvsd_orig = zeros(iter,1);
cvrms_orig = zeros(iter,1);
dfa_orig = zeros(iter,1);

% to hold variables after inducing sampling error
cvsd_err = zeros(iter,1);
cvrms_err = zeros(iter,1);
dfa_err = zeros(iter,1);

for i = 1:iter
    % create the noisy data, using alpha as exponent
    
    x0 = spatialPattern([N-1 1],alpha);
    x=zeros(N,1);
    x(1) = 0; % make first element 0 ms
    x(2:N) = x0; 
    x = (x - mean(x)) / std(x); % x now has mean = 0 and STD = 1
    
    % now rescale to target M and CV
    tarS = CV / 100 * M;
    x = x*tarS + M;
    
   % OK, now what happens when we induce error in sampling, due to device
   % sampling rate?
   sr = (1000 / Hz); % samples every xx ms
   
   x2 = sum_ser(x); % first, to be accurate, we need to recreate the actual time series; sum_ser is code by rje
   x3 = ceil(x2 / sr) * sr; % the new time stamp is rounded up to the next sampling point
   xx(1) = x(1);
   xx(2:N) = diff(x3); % this is now the resampled ITI series
   xx=xx(:);
   
   % use SMOOTH to see if we can recover the original seres
   % xxx = smooth(xx,0.1,'rloess'); % 0.1 is a 10% span; span must be odd
   xxx = smooth(xx,5);
    if iter == 1
        figure(100)
        subplot(3,1,1)
        plot(x,'b')
        hold on
        plot(xx,'r')
        plot(xxx,'g')
        hold off
        
        subplot(3,1,2)
        plot(x,xx,'.')
        
        subplot(3,1,3)
        plot(x,xxx,'.')
    end

    [a b c d] = rmssd(x);
    [a bb c dd] = rmssd(xx);
    
    cvrms_orig(i) = b;
    cvrms_err(i) = bb;
    cvsd_orig(i) = d;
    cvsd_err(i) = dd;
    
    % DFA
   [D alpha1] = DFA_main(x,4,60,20); % 4 to 60 is standard for gait analyis
   [D alpha2] = DFA_main(xx,4,60,20); % 4 to 60 is standard for gait analyis
   dfa_orig(i) = alpha1;
   dfa_err(i) = alpha2;

end

% relationship between CVSD and CVRMS

if iter > 1
    figure(300)
    subplot(1,3,1)
    plot(cvrms_orig,cvsd_orig,'.')
    xlabel('CV RMS'); ylabel('CV SD');

    subplot(1,3,2)
    plot(cvrms_orig,dfa_orig,'.')
    xlabel('CV RMS'); ylabel('DFA alpha');

    subplot(1,3,3)
    plot(cvsd_orig,dfa_orig,'.')
    xlabel('CV SD'); ylabel('DFA alpha');
    
    figure(301)
    subplot(1,3,1)
    plot(cvrms_orig,cvrms_err,'.');
    
    samp_errCV = cvrms_orig - cvrms_err; % or could do a bland altman plot
    samp_errSD = cvsd_orig - cvsd_err;
    
    subplot(1,3,2)
    hist(samp_errSD,20)
    xlabel('CVSD error')
    
    subplot(1,3,3)
    hist(samp_errCV,20)
    xlabel('CVRMS error')
    
end

if iter == 1
   output.cvsd_orig = cvsd_orig;
   output.cvsd_err = std(xx) / mean(xx) * 100;
   output.cvsd_sm = std(xxx) / mean(xxx) * 100;
   output.cvrms_orig = cvrms_orig;
   output.dfa_orig = dfa_orig;
   output.corr_orig_err = corr(x,xx);
   output.corr_orig_sm = corr(x,xxx);
   output.series = x;
else
   output.corr_cvrms_dfa = corr(cvrms_orig,dfa_orig);
end
