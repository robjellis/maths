function [output] = rms_sim(niter,win)

% base rates (per minute)
ratebpm = 50:25:150;
rates = 60000./ratebpm; % in ms
nrate = numel(rates);

% sampling (noise) rates (per second)
sampHz = [50 75 100 250 500];
noises = 1000./sampHz; % in ms
nnoise = numel(noises);

rmsCV = zeros(niter,nrate);

figure(100)

    for j = 1:nnoise
        for k = 1:nrate
            
            for i = 1:niter 
            nerr = 0 + (noises(j)-0).*rand(win,1); % uniform dist between 0 ms and noises(k) ms
            ev = zeros(win,1) + rates(k) + nerr;
            
            [rms rmscv] = rmssd(ev);
            
            rmsCV(i,k) = rmscv;
                     
            end  
        end
        
        % now get summary stats
        
        mnRMS = mean(rmsCV);
        plow = prctile(rmsCV,5);
        phigh = prctile(rmsCV,95);
        %errorbar(rates,mnRMS,plow,phigh); % plot with the 5 and 95 percentiles as error bars
        plot(ratebpm,phigh,'.-'); % just plot the UPPER percentile
        hold on
    end

hold off
xlabel('Target rate BPM')
ylabel('95th percentile CV_RMS')




output.mnRMS = mnRMS;