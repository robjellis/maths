function [seq] = geospace(a,b,n)

% function [seq] = geospace(a,b,n)
%
% create a geometric sequence from a to b in n steps, inclusive
%
% note: the "logspace.m" function (built in) does something different; this
% function actually returns "useable" values
%
% by rje

% make sure we go from small to large

    par = [a b];
    par = sort(par); % to make sure we have a < b
    a = par(1);
    b = par(2);

mult = (b / a)^(1/(n-1));

seq = a*mult.^(0:(n-1));

