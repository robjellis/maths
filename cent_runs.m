function [output] = cent_runs(X)

% mean center each run within the matrix X

nscans = size(X,1);
nrois = size(X,2);

fprintf(strcat(['\n There are ', num2str(nscans), ' scans (rows).']));
nruns = input('\n Number of runs?: ');
spr = input(' Number of scans per run?: ');

output = nan(size(X));

inc = 1; % start with row 1

for r = 1:nruns
    x = X(inc:(inc+spr-1),:);
    mx = mean(x);
    mx = mean(mx);
    
    for j = 1:nrois
        output(inc:(inc+spr-1),j) = x(:,j) - mx;
    end
    inc = inc + spr;
    
end


