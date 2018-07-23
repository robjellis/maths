% mct
%
% Monte Carlo simulations for pseudo T distributions to see the range of
% skewness or T_crit values that appear.
%
% Copyright 2010 Robert J Ellis

fprintf('\n Monte Carlo simulations of Student T-distributions. \n\n')

df = input(' Enter the desired d.f.: ');
N = input(' Enter the number of obs. in the distribution: ');


    alp = input(' Enter the voxel-level alpha to determine: ');
    percval = 100*(1-alp);
    t_crit = tinv(1-alp,df);

perm = input(' Enter the number of M.C. permutations: ');

fprintf(['\n Running simulations ... \n']);

skew = zeros(perm,1);
perc = zeros(perm,1);
Avox = zeros(perm,1);
Mt = zeros(perm,1);
Ap = zeros(perm,1);

% makes use of progressbar.m and gui_active.m by Ohad Gal, obtained from MATLAB central
% at http://www.mathworks.com/matlabcentral/fileexchange/3607-progressbar

gui_active(1);      % will add an abort button
h = progressbar( [],0,'Running M.C. simulations ...' );
    
for i = 1:perm
    tdist = trnd(df-1,N,1);

        skew(i) = skewness(tdist);

        perc(i) = prctile(tdist,percval);
        Avox(i) = sum(tdist>t_crit);
        Ap(i) = Avox(i)/N;
        Mt(i) = mean(tdist);
   
    
    h = progressbar( h,1/perm );                 % the progress bar 
    
    if ~gui_active
        break;
    end
end

progressbar( h,-1 );



    skval = num2str(prctile(skew,95));
    simval = num2str(prctile(perc,95));
    simval2 = num2str(prctile(Avox,95));
    simval3 = num2str(prctile(Mt,95));
    simval4 = num2str(prctile(Ap,95));
        
    fprintf(['\n\n M.C. simulation complete. \n 95th percentile of "T_crit" distribution (voxel alpha = ' num2str(alp) ') is ' simval '.\n']);
    fprintf([' 95th percentile of "voxels > T_crit" distribution is ' simval2 '.\n']);
    fprintf([' 95th percentile of "mean T-value" distribution is ' simval3 '.\n']);
    fprintf([' 95th percentile of "P(i+)" distribution is ' simval4 '.\n']);    
    fprintf([' 95th percentile of "skewness" distribution is ' skval '.\n']);
    fprintf([' (The theoretical T_crit at alpha = ' num2str(alp) ' is ' num2str(t_crit) '.) \n\n']);



    