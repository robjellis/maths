function output = propsup(s1,s2)

% proportion of superiority, manually calculated

n1 = numel(s1);
n2 = numel(s2);

props = nan(n1,1);

for i = 1:n1
   props(i) = sum(s1(i) > s2) / n2;   
end

% get the PS
output.PS = mean(props);
