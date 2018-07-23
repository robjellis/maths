function ex_analysis(ex_inds)

% input ex_inds and get some plots as output

numf = numel(ex_inds);

x = 10:250; % total possible footsteps

all_y = zeros(numf,numel(x));

for f = 1:numf
    inds = ex_inds{f};
    %numi = numel(inds);
    %inds_norm = inds / numi * 100; % percentage
    
    y = histc(inds,x);
    all_y(f,:) = y';
    
end

all_y = sum(all_y) / numf * 100;

figure(30)
plot(x,all_y)
ylabel('Percentage of total files')
xlabel('Ordinal step')
title('Positions of deleted steps (turns)')
