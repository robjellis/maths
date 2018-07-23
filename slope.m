function output = slope(x1,y1,x2,y2,x3)

% calculate y3 based upon the two other coordinates


%to find a target x (beta) value based on two x,y coordinate pairs and a target y (success probability) value 

slope = (y2-y1)/(x2-x1);
b = y1 - slope*x1;
y3 = slope*x3 + b;

figure(10)
plot(x1,y1,'b.')
hold on
plot(x2,y2,'b.')
plot(x3,y3,'r.')
hold off

output.slope = slope;
output.intercept = b;
output.x3 = x3;
output.y3 = y3;



    