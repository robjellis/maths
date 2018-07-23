% mct
%
% Monte Carlo simulations for pseudo T distributions to see the range of
% skewness or T_crit values that appear.
%
% Copyright 2010 Robert J Ellis

fprintf('\n Monte Carlo simulations of Student T-distributions. \n\n')

df = input(' Enter the desired d.f.: ');
N = input(' Enter the number of observations (i.e., voxels) in the distribution: ');
perm = input(' Enter the number of M.C. permutations: ');
percdist = input(' Enter the percentile value of the M.C. distribution: ');

fprintf(['\n Running simulations ... \n']);

skew = zeros(perm,1);
perc = zeros(perm,1);
Avox = zeros(perm,1);
Mt = zeros(perm,1);
Ap = zeros(perm,1);
pcorr = zeros(perm,1);

% makes use of progressbar.m and gui_active.m by Ohad Gal, obtained from MATLAB central
% at http://www.mathworks.com/matlabcentral/fileexchange/3607-progressbar

gui_active(1);      % will add an abort button
h = progressbar( [],0,'Running M.C. simulations ...' );
    
for i = 1:perm
    tdist = trnd(df,N,1);
    tdist2 = trnd(df,N,1);

        skew(i) = skewness(tdist);
      
        pcorr(i) = corr(tdist,tdist2);
   
    
    h = progressbar( h,1/perm );                 % the progress bar 
    
    if ~gui_active
        break;
    end
end

progressbar( h,-1 );



    skval = num2str(prctile(skew,percdist));

    corrval = num2str(prctile(pcorr,percdist));
    percdistv = num2str(percdist);
    fprintf('\n\n    ----------------------------  ');
    fprintf('\n\n    M.C. simulation complete. \n\n');  
    fprintf(['    ' percdistv ' percentile of "skewness" distribution is ' skval '.\n']);
    fprintf(['    ' percdistv ' percentile of "r-value" distribution is ' corrval '.\n\n']);
    fprintf('\n    ---------------------------- \n '); 



    