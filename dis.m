% 3-D Euclidian distance (for brain coordinates)
%
% Prompts user to enter two points [x1 y1 z1] and [x2 y2 z2] and calculates
% distance.
%
% rje 2012

function [d] = dis


a = input(' Enter 1st point [x y z]: ');
b = input(' Enter 2nd point [x y z]: ');

x1 = a(1); y1 = a(2); z1 = a(3);
x2 = b(1); y2 = b(2); z2 = b(3);

d = sqrt((x2 - x1)^2 + (y2 - y1)^2 + (z2 -z1)^2);