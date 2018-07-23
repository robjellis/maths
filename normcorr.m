function [r NC] = normcorr(X,Y)

% regular correlations

r = corr(X,Y);

NC = sum(X.*Y) / (sqrt(sum(X.^2))*sqrt(sum(Y.^2)));