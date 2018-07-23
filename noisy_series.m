function [output x xx] = noisy_series(M,CV,beta,iter,N,Hz)

%
% function [output] = noisy_series(M,CV,beta,iter,N,Hz)
%
% create a noisy series, and see how outcome measures are affected by
% changes in "sampling rate"
%
% [output] = noisy_series(iter,N,M,CV,beta,Hz)
%
% outside the iteration loop
%    M     = true seq mean (in ms) [x-axis; dim a]
%    CV    = true CV (percentage) [sep lines; dim b]
%    beta  = noise color [sep plots; dim c]
%            0 is random white noise  
%           -1 is pink noise
%           -2 is Brownian noise 
%
% then we do the iteration, creating a new sequence with the above
% properties
%
% inside each iteration, then we have:
%    N     = length of sequence [dim d]
%    Hz    = "observed" sampling rate, in ms (e.g., 100 is coarser than 1000) [dim e]
%
% Figure has three plots:
%   1. Blue is original ITI series; red is the series with noise
%   2. Correlation of Original ITIs with Noisy ITIs
%   3. Poincare plot of the Original series
%
% rje | version = 2013.03.01
%

% calls function from:
% http://www.mathworks.com/matlabcentral/fileexchange/5091-generate-spatial-data/content/spatialPattern.m
%

%% outcome measures


num_d = numel(N);
num_e = numel(Hz);

max_N = max(N); % we always start with the most number of elemens

% these are all performed, for each sequence, by td_calcs (which outputs in
% rows)


for a = 1:numel(M)
    for b = 1:numel(CV)
        for c = 1:numel(beta)
            
            progressbar % for iteration loop
            for i = 1:iter
                %% create the noisy data, using beta as exponent
                progressbar(i/iter)
                if beta == 0
                    % just use true randn sequence
                    x0 = randn(max_N-1,1);
                else 
                    x0 = spatialPattern([max_N-1 1],beta); % http://www.mathworks.com/matlabcentral/fileexchange/5091
                end
                x = (x0 - mean(x0)) / std(x0); % x now has mean = 0 and STD = 1

                % now rescale to target M and CV
                tarS = CV(b) / 100 * M(a);
                x = x*tarS + M(a);

                % finally, get the "actual" time series
                xt(1) = 0; % make the first timestamp at 0 ms
                xt(2:max_N) = cumsum(x);

                %% induce sampling "error" 
                % what happens when we induce error in sampling, due to device sampling rate?

               sr = (1000 / Hz); % samples every xx ms

               xte = ceil(xt / sr) * sr; % the new time stamp is rounded *up* to the next sampling point

               xx = diff(xte); % this is now the resampled ITI series, which is compared to original x series
               xx=xx(:);

               %mad(xx,1)
               %std(xx)

               % use SMOOTH to see if we can recover the original seres
               % xxx = smooth(xx,0.1,'rloess'); % 0.1 is a 10% span; span must be odd
               % xxx = smooth(xx,5); % moving average of 5 events

                if iter == 1

                    % just to make an easier plot
                    figure(99)
                    plot(x,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0],'MarkerSize',15,'Marker','.','LineWidth',1,'Color',[0 0 1])
                    xlabel('Ordinal Inter-step Intervals')
                    ylabel('Inter-step Interval (ms)')

                    figure(100)
                    subplot(3,1,[1 2])
                    plot(xx,'r','LineWidth',1) % plot this behind
                    hold on
                    plot(x,'b','LineWidth',1) % original
                    hold off

                    subplot(3,1,3)
                    plot(x,xx,'.')
                    xlabel('Original ITIs'); ylabel('ITIs + sampling noise');

                    %subplot(4,1,4)
                    %plot(x(1:end-1),x(2:end),'.')
                    %xlabel('1 to max_N - 1'); ylabel('2 to max_N');

                end

                % use td_calcs to do all the time-domain calculations 
                td_o = td_calcs( x,[],50:5:95,0); 
                td_e = td_calcs(xx,[],50:5:95,0);



            end % iter loop
        end % c loop
    end % b loop
end % a loop


%% plots

if iter > 50
    % relationship between CVSD_e and CVRMS_e
    plot_eda(cvsd_e,cvrms_e,9,50,'CVSD_e','CVRMS_e')
      
end


% get percentiles
ll = 5;
mm = 50;
hh = 95;

if iter == 1
   output.sd_o = sd_o;
   output.sd_e  = sd_e;
   output.cvsd_o = cvsd_o;
   output.cvsd_e = cvsd_e;
   %output.cvsd_sm = std(xxx) / mean(xxx) * 100;
   output.rms_o  = rms_o;
   output.rms_e   = rms_e;
   output.cvrms_o = cvrms_o;
   output.cvrms_e = cvrms_e;
   output.dfa_o = dfa_o;
   output.dfa_e = dfa_e;
   output.corr_o_e = corr(x,xx);
   %output.corr_o_sm = corr(x,xxx);
   
   output.poin_rho_o = poin_rho_o;
   output.poin_rho_e = poin_rho_e;
   output.orig_series = x;
   output.noisy_series = xx;
else
    output.target_M = M;
    output.target_CV = CV;
    output.target_alpha = beta;
    %output.cvsd_e = cvsd_e;
    output.cvsd_e_minus = prctile(cvsd_e,mm) - prctile(cvsd_e,ll);
    output.cvsd_e_med = prctile(cvsd_e,mm);
    output.cvsd_e_plus = prctile(cvsd_e,hh) - prctile(cvsd_e,mm);
    output.cvrms_o = cvrms_o;
    output.cvrms_e = cvrms_e;
    output.dfa_o = dfa_o;
    output.dfa_e = dfa_e;
    output.dfa_o_minus = prctile(dfa_o,mm) - prctile(dfa_o,ll);
    output.dfa_o_med = prctile(dfa_o,mm);
    output.dfa_o_plus = prctile(dfa_o,hh) - prctile(dfa_o,mm);
    output.dfa_e_minus = prctile(dfa_e,mm) - prctile(dfa_e,ll);
    output.dfa_e_med = prctile(dfa_e,mm);
    output.dfa_e_plus = prctile(dfa_e,hh) - prctile(dfa_e,mm);
    output.corr_dfa_o_e = corr(dfa_o,dfa_e);
end
