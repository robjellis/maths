function [output] = mc_prc_meas(N,M,CV,outperc,iter)

% [output] = mc_prc_meas(N,M,CV,outlier,iter)
%
% Monte Carlo simulations of how the RJE-defined time domain percentile
% measures (in td_calcs.m) perform relative to standard measures (SD and
% RMSSD) under pure random noise conditions with a fixed percentage of
% outliers
%
% "outperc = X" will randomly (linspace) replace X% of the data with high value outliers 
%
%
% RJE | 2013.03.06

prc = [5 : 5 : 95]; % percentiles to evaluate

nprc = numel(prc);

% variables to save - will either be raw or CV format
  sd = nan(iter,1);
 rms = nan(iter,1);
padm = nan(iter,numel(prc));
pasd = nan(iter,numel(prc));

for n = 1:numel(N)
    
    for c = 1:numel(CV)
        
        for i = 1:iter
            
            if CV == 0
                x = randn(N(n),1)+M; % SD = 1;
            elseif CV > 0
                % new vector of ITIs, based upon a 1000-ms period
                % we can't get a valid CV if M = 0!
                x = randn(N(n),1);
                tar_sd = (CV(c) * M) / 100 ;
                %mult = tar_sd / std(x);  % we can't get this to be exact, or else we can't do a correlation!

                x = x*tar_sd; % note: the observed SD will vary from run to run

                x = x - mean(x) + M; % now will have the correct mean and CV
            
            end
            
            sdx = std(x);
            
            if outperc > 0
                nout = round(outperc/100 * numel(x));
                xpos = round(linspace(1,numel(x),nout+2)); % the positions
                xpos(1) = [];
                xpos(end) = []; % don't use the first and last positions

                x(xpos) = mean(x) + 10 * sdx; % high value outlier
                
                %figure(11)
                %plot(x)
  
            end
            
            res = td_calcs(x,prc);          

            sd(i)     = res.cvsd;
            rms(i)    = res.cvrms;
            padm(i,:) = res.padm_raw;
            pasd(i,:) = res.pasd_raw;

        end
    end
end

%% summary measures
% Use Bland-Altman SD error values as the summary statistic, using rje code
% "ba_calc.m"

% loops
sd_padm  = nan(2,nprc); % row1 = mean; row2 = 2SD
rms_pasd = nan(2,nprc);
adm_spearman   = nan(1,nprc);
asd_spearman   = nan(1,nprc);
for p = 1:nprc
     [res] = ba_calc(sd,padm(:,p),0); % because B-A plot is [meas1 - meas 2] for y-axis, positive values will mean that sd > PADM
     
     %figure(30)
     %plot(sd,padm(:,p),'.')

     sd_padm(1,p) = res.BA_mean;
     sd_padm(2,p) = res.BA_2sd;
     adm_spearman(1,p) = res.Pearson; % Spearman is too forgiving in this instance
     
     [res] = ba_calc(rms,pasd(:,p),0);
     rms_pasd(1,p) = res.BA_mean;
     rms_pasd(2,p) = res.BA_2sd;
     asd_spearman(1,p) = res.Pearson;
end

output.percentiles = prc;
output.sd_padm  = sd_padm;
output.rms_pasd = rms_pasd;


    figure(10)
    
    subplot(1,3,1) % Spearman correlation
        plot(prc,adm_spearman,'b','MarkerSize',10,'Marker','.','DisplayName','SD vs. PADM')
    hold on
        plot(prc,asd_spearman,'r','MarkerSize',10,'Marker','.','DisplayName','RMS vs. PASD')
    hold off
    
    xlabel('Percentile')
    ylabel('Pearson r-value')
    legend('show','Location','NorthWest');
    
    subplot(1,3,2) % plot B-A mean
        plot(prc,output.sd_padm(1,:),'b','MarkerSize',10,'Marker','.','DisplayName','CVSD vs. PADM')
    hold on
        plot(prc,output.rms_pasd(1,:),'r','MarkerSize',10,'Marker','.','DisplayName','RMS vs. PASD')
    hold off
    
    xlabel('Percentile')
    ylabel('Bland-Altman Mean')
    legend('show','Location','NorthWest');
    
    subplot(1,3,3) % plot B-A 2SD
        plot(prc,output.sd_padm(2,:),'b','MarkerSize',10,'Marker','.','DisplayName','SD vs. PADM')
    hold on
        plot(prc,output.rms_pasd(2,:),'r','MarkerSize',10,'Marker','.','DisplayName','RMS vs. PASD')
    hold off
    
    xlabel('Percentile')
    ylabel('Bland-Altman 2SD')
    legend('show','Location','NorthWest');
