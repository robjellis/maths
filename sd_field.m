% error "field" statistic

% dummy variables for now
a = hhh(:,1);
b = hhh(:,2);
plot(a,b,'.')
step = .01;
div = 1/step;

% round along x dimension to 2 decimal places
a = round(a(:,1)*div)/div;

amin = min(a(:,1));
amax = max(a(:,1));
xrnd = xmin:step:xmax;
alength = length(xrnd);

% now find the std dev of y values for each rounded x value

for i=1:alength
    res(i,1) = xrnd(i);
    mn(i) = nanmean(b(a == xrnd(i)));
    sd(i) = nanstd(b(a == xrnd(i)));
    res(i,2) = mn(i) - sd(i);
    res(i,3) = mn(i) + sd(i);
    
    
end

res;
hold on
plot(res(:,1),res(:,2),'r')
plot(res(:,1),res(:,3),'g')
hold off

    
res
res = []

