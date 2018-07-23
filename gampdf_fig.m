function output = gampdf_fig(k,theta,M,CV)

% this to generate Y-values (Gamma PDF values) for a given k value, target
% mean, and target CV, so it makes it easier to compare the shapes of gamma
% with respect to gait or HR data.
%
% Will loop through multiple levels of CV and add them to a plot

%% output vars to save
tarn = [10^2 10^3 10^4 10^5];
tarn = tarn(:);

% now get the p-value to check
tarp = sqrt(1./tarn);

n_vs_rat = nan(numel(tarn),numel(CV+1));
n_vs_rat(:,1) = tarn; % so we have "headers"

%% do it
% get the xvalues
pmin = 10 ^(-10);

% get xmin and xmax
xmin = gaminv(pmin,k,1); 
xmax = gaminv(1-pmin,k,1);

x = linspace(xmin,xmax,10^(4)); 
x = x(:);

% now get the PDF values
pdf = gampdf(x,k,M);

% what is the "estimated" mean, based on a simple calulation from the PDF?

est_mean = sum(x.*pdf) / sum(pdf);

%% add the CV loop here

for c = 1:numel(CV)
    % now we have to rescale the x-values; borrows logic from rje smooth_series
    theoM = k*theta;
    theoS = sqrt(k*theta^2);
    tarS = (CV(c) * M) / 100;

    xmod = (x - theoM * 1000*M)/1000; % now it will have target M, but not target SD

    xmod = (xmod / (theoS/1000)) * tarS; % theoS/1000 is critical because we divided the whole series by 1000

    % what is the *actual* mean of this? Well, we can't simply take the mean of
    % the xmod values, since the PDF is not uniform! So, we have to approximate
    % this using numerical integration:
    % http://mathworld.wolfram.com/NumericalIntegration.html; i.e., we take the
    % weighted mean of all the rectangles; the more rectangles we have (and the
    % better the range they cover, the better)

    xmod_mean_est = sum(xmod.*pdf) / sum(pdf); % the more complete the x series, the better this value will be; RJE confirmed this

    % now the mean will be wrong, so we adjust the mean one more time
    xmod = xmod - xmod_mean_est + M;

    mean_rs = sum(xmod.*pdf) / sum(pdf); % should always be 1.0

    % also may be helpful to rescale pdf values such that the largest value =
    % 1.0

    pdf_rs = pdf / max(pdf);

    %% now we create the CDF, and get the desired P-values

    cdf = cumsum(pdf) / max(cumsum(pdf));

    for n = 1:numel(tarn)
        lowV =  max(xmod(cdf <= tarp(n)));
        highV = min(xmod(cdf >= (1 - tarp(n))));

        n_vs_rat(n,c+1) = highV / lowV; % this is the largest ratio
    end
    
    
    figure(11)
    plot(xmod,cdf)
    hold on





    %% figures

    figure(10)
    plot(xmod,pdf)
    hold on

end
figure(10)
hold off

figure(11)
hold off

%% output 

output.x = x;
output.x_rs = xmod;
output.pdf = pdf;
output.pdf_rs = pdf_rs;
output.orig_est_mean = est_mean;
output.rs_est_mean = mean_rs;
output.n_vs_rat = n_vs_rat;

