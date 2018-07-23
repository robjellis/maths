function [PSD,F]=welch_rje(t,y,window,noverlap,nfft,fs)

%Prepare y
t2 = t(1):1/fs:t(length(t));%time values for interp.
y=interp1(t,y,t2','spline')'; %cubic spline interpolation
y=y-mean(y); %remove mean

%Calculate Welch PSD using hamming windowing    
[PSD,F] = pwelch(y,window,noverlap,(nfft*2)-1,fs,'onesided');

%% plot
figure(200)
plot(F,PSD)