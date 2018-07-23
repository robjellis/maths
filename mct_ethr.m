% mct
%
% Monte Carlo simulations for pseudo T distributions to see the range of
% skewness or T_crit values that appear.
%
% Note: dfd and VB (voxels in brain) are obtained from spm_tmax
%
% Copyright 2010 Robert J Ellis

perm = input(' Enter the number of M.C. permutations: ');
percdist = input(' Enter the percentile value of the M.C. distribution: ');
plotfigs = input(' Display M.C. results? [y/n]: ','s');


skew = zeros(perm,1);
perc = zeros(perm,1);
pcorr = zeros(perm,1);

histstep = -6:.25:6;
counts = zeros(perm,numel(histstep));

% makes use of progressbar.m and gui_active.m by Ohad Gal, obtained from MATLAB central
% at http://www.mathworks.com/matlabcentral/fileexchange/3607-progressbar

gui_active(1);      % will add an abort button
h = progressbar( [],0,'Running M.C. simulations ...' );


    
for i = 1:perm
    tdist = trnd(dfd,VB,1);
    tdist2 = trnd(dfd,VB,1);
    counts(i,1:numel(histstep)) = histc(tdist,histstep)';
    
        skew(i) = skewness(tdist);
      
        pcorr(i) = corr(tdist,tdist2);
   
    
    h = progressbar( h,1/perm );                 % the progress bar 
    
    if ~gui_active
        break;
    end
end

progressbar( h,-1 );

% === plot simulation results

if plotfigs == 'y'
    
figure
for i = 1:perm
plot(histstep,counts(i,:));
hold on
end

figure
hist(skew,20)

figure
hist(pcorr,20)
end

    skval = prctile(skew,percdist);
    corrval = prctile(pcorr,percdist);
    percdistv = num2str(percdist);
    
    fprintf('\n    ----------------------------  ');
    fprintf('\n    M.C. simulation complete. \n\n');  
    fprintf(['    ' percdistv ' percentile of "skewness" distribution is ' num2str(skval) '.\n']);
    fprintf(['    ' percdistv ' percentile of "r-value" distribution is ' num2str(corrval) '.']);
    fprintf('\n    ---------------------------- \n '); 



    