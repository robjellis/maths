function out = mom_beta(x,s,scale,plot_it)

% method of moments for beta distribution
% input observed data mean(x) and sample SD (s); x values must be 0 to 1 scaled
% the corresponding estimated alpha and beta params will be returned
%
% RJE | 21 Feb 2018
%
% x and s can be vectors; NaNs are acceptable and will return as NaN
%
% example references:
% - http://www.itl.nist.gov/div898/handbook/eda/section3/eda366h.htm 
% - http://www.real-statistics.com/distribution-fitting/method-of-moments/method-of-moments-beta-distribution/

if numel(x) == 1
    plot_it = 1;
end

% values must be on a 0 to 1 scale; rescale if not

if scale == 1
    
elseif scale == 5
    % 1 to 5 scale linear transformed to 0 to 1 scale
    x = 0.25 * x - 0.25;
    s = s * 0.25;
else
    error(' Values must be on a 0 to 1 scale')
end

% square s for easier typing
s2 = power(s,2);

a_hat = x .* ((x .*(1-x)./s2) - 1);

b_hat = a_hat .* (1 - x) ./ x;

out.a_hat	= a_hat;
out.b_hat	= b_hat;

if plot_it == 1 
    figure(40)
    xval = 0:.001:1;
    yval = betapdf(xval,a_hat,b_hat);
    plot(xval,yval,'r')
end


