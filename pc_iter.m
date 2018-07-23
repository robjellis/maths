% evaulting partial correlations in a regression
% default: two regressors (x1 and x2) and one dependent (y)
% robert j ellis feb 2011

function pc_iter(x1,x2,yvals)

iter = 10000;
lenx1 = numel(x1);
lenx2 = numel(x2);


x1 = x1(:); x2 = x2(:);  % just get everything in columns

pc_x1_y = zeros(iter,1);
pc_x2_y = zeros(iter,1);

if lenx1 ~= lenx2
   fprintf(' Your regressor vectors do not have the same length! ');
   return
end

options = input(' Options: \n [1] random y values, unyoked; [2] random y values, yoked with Age; \n [3] random y values, yoked with Training; \n [4] permutations of y-vals: ');
%matchit = input(' Match input to [1] Age or [2] Training?: ');

    
for i = 1:iter
    
    if options == 1 || options == 2 || options == 3
    y = randn(lenx1,1);  % new random vector on every trial
    end
    if options == 2 || options == 3
    y = sort(y);
    end

    if options == 2
       x1 = sort(x1);
       y = (y .* x1); 
    elseif options == 3
       x2 = sort(x2);
       y = (y .* x2);  
    end
    if options == 4
        y = zeros(lenx1,1);
    sortme = randperm(lenx1);
    for j = 1:lenx1
       y(j) = yvals(sortme(j)); 
    end
      y;
    end
    
    x1_y = [x1 y];
    x2_y = [x2 y];
    
    temp1 = partialcorr(x1_y,x2); % partial corr of x1 and y, controlling for x2
    pc_x1_y(i) = temp1(1,2);
    
    temp2 = partialcorr(x2_y,x1); % partial corr of x2 and y, controlling for x1
    pc_x2_y(i) = temp2(1,2);
        
end

% make the figure

x1vals = linspace(min(pc_x1_y),max(pc_x1_y),30);
x2vals = linspace(min(pc_x2_y),max(pc_x2_y),30);

x1hist = histc(pc_x1_y,x1vals);
x2hist = histc(pc_x2_y,x2vals);
figure
plot(x1vals,x1hist,'b');
hold on
plot(x2vals,x2hist,'r');
hold off

mn_pc_x1_y = mean(pc_x1_y) 
std_pc_x1_y = std(pc_x1_y)
mn_pc_x2_y = mean(pc_x2_y) 
std_pc_x2_y = std(pc_x2_y)

corr(pc_x1_y,pc_x2_y)
figure
plot(pc_x1_y,pc_x2_y, '.')
