function output = smooth_series(DTYPE,DPARAM,SEED,M,CV,SMTYPE,SM,Nin,Nout,PLOT_IT)

% output = smooth_series(DTYPE,DPARAM,SEED,M,CV,SMTYPE,SM,Nin,Nout,PLOT_IT)

% code to make a *single* smooth series with the above input parameters
%
% note: input mean in *seconds*, not ms, to facilitate plotting
%
% SEED is a non-negative integer; can be [] to ignore
% Nin is *internal* length of the full sequence to generate
% Nout is the *external* length of the output sequence 


%% color or black and white
use_color = 0;

if use_color == 0
    gtcol = [0.502 0.502 0.502]; % gray
    dscol = [0 0 0]; % black
elseif use_color == 1
    gtcol = 'b';
    dscol = 'r';
end

%% first we generate the uniform dist so we can then smooth it

% commmon to all methods
if isempty(SEED)
    % just regular rand  
else 
    rand('twister',SEED) % will be a one-time-use only 
end

if strcmp(SMTYPE,'ma')
    % just add the width of the window on either side
    XTRA = SM;
    po = rand(Nin+2*XTRA,1);
    
    ps = smooth(po,SM,'moving');
    % Note from Matlab: We don't need to specify "x". 
    % When X is given and X is not uniformly distributed, the default method is 'lowess'.  
    % The 'moving' method is not recommended.
        
elseif strcmp(SMTYPE,'sp')
    XTRA = 50;
    po = rand(Nin+2*XTRA,1);
    
    % still have to write this
    
elseif strcmp(SMTYPE,'wp') % by rje: "weighted priors"
    % number of events in the prior
    WIN = 20; % rje: 20 is best for the current (2013-08-21) method of smoothing the CDF and then getting the distribution values, 
              % in order to yield a nice spread for DFA
    
    XTRA = WIN; % we only look back at this many steps
    po = rand(Nin+2*XTRA,1);
       
    w = SM; % the weight

    for j = WIN:-1:1
        if j == WIN
           weights(j) = 1;
        else
           weights(j) = weights(j+1)*w;
        end
    end
    %weights = w.^[0:(WIN-1)];
    weights = weights(:); % a column to match data
    sum_weights = sum(weights);

    % note: we will ignore the first WIN steps (just copy values for convenience

    ps = nan(numel(po),1);
    for i = 1:numel(po)
        if i < WIN
            ps(i) = po(i);
        else
            y_sec = po((i-WIN+1) : i);
            ps(i) = sum(y_sec.*weights) / sum_weights;
        end
    end
end % smoothing type

% ************************
% now we proceed similarly for all methods

% now we just take the middle smoothed portion; these are now treated as the the y-values on the CDF
poc = po(XTRA+1:XTRA+Nin); % still on the Nin values
psc = ps(XTRA+1:XTRA+Nin);

% first, let's just replace extreme values in poc
poc(poc > .995) = .995;
poc(poc < .005) = .005;

% rescale all values back to target uniform distribution
tarS = sqrt(1/12);
tarM = 0.5;

psc = psc - mean(psc); % will have mean of 0

% now we have to figure out the multiplier so that we we restore the min or
% max value from the original series (minus 0.5), whichever is smaller. In
% this way, we get more variety in the values that we see; also, the longer
% the data set, the closer to the theoretical min and max we will be!

% revision: actually, we start to run into problems when we get values that
% approch .001 or .999, especially when CV is large. But, that is actually
% more accurate! We can also now run our outlier detection method on the
% simulated data to see how much data we cut out!

mult =  min((max(poc) - 0.5) / max(psc), (min(poc) - 0.5) / min(psc)); % original will always be larger than smoothed

% now multiply by this value; mean will still be at 0
psc = psc * mult;

%shift up to the target mean of 0.5
psc = psc - mean(psc) + 0.5;

% old - don't need to do this!
% psc(psc >= max(poc)) = max(poc);
% psc(psc <= min(poc)) = min(poc);

% note: the target SD will not be met, but we really can't control that and
% still have the range we want!

% mean(psc)
% std(psc)
% 
% 



%% now we get the target distribution values using inverse CDF functions
% we also rescale to the target M and CV

% get the target SD - works same for all methods
SD = (M * CV)/100;

if strcmp(DTYPE,'u')
    x = psc;
    % rescale to target
    x = (x - mean(x)) / std(x); % the center section now has mean = 0 and STD = 1

    x = x*SD + M;
    
elseif strcmp(DTYPE,'n')
    % there is no distribution parameter at this point, so DPARAM should = [];
    x = norminv(psc,0,1);
    
    % rescale to target
    x = (x - mean(x)) / std(x); % the center section now has mean = 0 and STD = 1

    x = x*SD + M;
    
elseif strcmp(DTYPE,'p')
    % Poisson distribution; has a minimum of zero
    lambda = DPARAM;
    x = poissinv(psc,lambda);    

    % still need to rescale
    
elseif strcmp(DTYPE,'g')
    % gamma distribution; has a minimum of zero
    shape = DPARAM;
    scale = 1; % keep this as 1 so that mean ~ shape and std = sqrt(shape)
    x = gaminv(psc,shape,scale);
    
    % 1. add the mean in *ms*
    x = x - mean(x) + 1000*M; % y will now have a precise mean = M*1000
    
    % 2. get back into seconds
    x = x / 1000;

    % 3. rescale to target SD; again, just use center section    
    x = (x / std(x)) * SD;
    
    % 4. now subtract the mean one more time
    x = x - mean(x) + M;
    
end

%% OK, a final recentering for Nout ...

xx = x(1:Nout);
xx = (xx - mean(xx)) / std(xx); % now has mean = 0 and STD = 1
xx = xx*SD + M;
    
if PLOT_IT == 1    
    
    %figure(20)
    %autocorr(psc)

    figure(100)
    subplot(3,1,1)
    pind_all = 1:numel(po);
    pind_cen = pind_all(XTRA+1:XTRA+Nin);

    plot(pind_all,po,'k','LineWidth',1,'Marker','.')
    hold on
    plot(pind_cen,poc,'r','LineWidth',1,'Marker','.')
    plot(pind_cen,psc,'b','LineWidth',1,'Marker','.')
    axis([1 numel(po) 0 1])
    xlabel('Index')
    ylabel('CDF y-value')
    title('(From smooth_series.m)','Interpreter','none')

    hold off

    subplot(3,1,2)
    % just the smoothed series

    plot(pind_cen,psc,'b','LineWidth',1,'Marker','.')
    xlabel('Index')
    ylabel('CDF y-value')
    axis([1 numel(po) 0 1])
    
    subplot(3,1,3)
    % real values
    plot(pind_cen,x,'b','LineWidth',1,'Marker','.')
    hold on
    xlabel('Index')
    ylabel('Distribution x-values')
    %plot(pind_cen(1:Nout),x(1:Nout),'r','LineWidth',1,'Marker','.')
    plot(pind_cen(1:Nout),xx,'r','LineWidth',1,'Marker','.')
    axis([1 numel(po) min(x) max(x)])
    hold off
    

end

% just output the center sections
output.cdf_orig_cent = poc; % original CDF center values
output.cdf_smth_cent = psc; % smoothed CDF center values
output.val_smth_Nin = x; % distribution values
output.val_smth_mean = mean(x);
output.val_smth_std = std(x);
output.val_smth_Nout = xx;

% let's output the correct time series, with a 0 in position 1, so we can
% re-create the exact IEI series during td_calcs_concat

yy = cumsum(xx);
yy(2:end+1) = yy;
yy(1) = 0;

output.val_smth_Nout_ts = yy;

    