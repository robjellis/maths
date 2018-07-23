function output = kde_sim(N,DTYPE,div,iter)

% testing the size of the mesh for kde.m

peaks = nan(iter,2);

for i = 1:iter
    if strcmp(DTYPE,'n')
        x = randn(N,1);
    elseif strcmp(DTYPE,'u')
        x = rand(N,1); % mean = 0.5
    end
    
    % default
    [bandwidth1, density1, xmesh1] = kde(x,N);
    
    % sparser
    [bandwidth2, density2, xmesh2] = kde(x,N/div);
    
    % get the peak
    peak1 = xmesh1(density1 == max(density1));
    peak2 = xmesh2(density2 == max(density2));
    
    peaks(i,1) = peak1; 
    peaks(i,2) = peak2;
    
end

plot_eda(peaks(:,1),peaks(:,2),5)
