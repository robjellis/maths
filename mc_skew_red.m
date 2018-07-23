% Monte Carlo tests of skewness for T-distributions
% 
% the underlying distribution is important (e.g., T-distribution), as
% adding high value voxels will make more of an impact when df is higher!

N = input(' Number of voxels in the brain: ');
df = input (' Enter the d.f.: ');
iter = input(' Number of Monte Carlo iterations: ');
prt = input(' Percentile to evalute (1 to 99): ');
makecorr = input(' Make correlation plots [y/n]? ','s');

pvals = [.0001 .00025 .0005 .001 .0025 .005];

% signal sizes
sig = ceil(N .* pvals);       % sizes will go from smallest to largest

% signal strength means
strength = abs(tinv(pvals,df));  % t-values will go from largest to smallest

gui_active(1);                                    % will add an abort button
h = progressbar( [],0,'Working ...' );          % progress bar


for k = 1:length(pvals)

    h = progressbar( h,(1/length(strength)));                 
    
    if ~gui_active
        break;
    end  
    
    
for i = 1:iter
    
  vol1 = trnd(df,N,1);
  vol1(1:sig(j)) = randn(sig(j),1) + strength(k);
  
  vol2 = trnd(df,N,1);
  vol2(1:sig(k)) = randn(sig(k),1) + strength(k);
  
  vol_sk(i, k) = skewness(vol1);
  vol_corr(i, k) = corr(vol1,vol2);

if makecorr == 'y' && i == iter && k == length(sig)
   figure
   plot(vol1,vol2,'.');
   xlabel('Values from vol1'); 
   ylabel('Values from vol2');
   title([num2str(sig(k)) ' voxels with M = ' num2str(strength(k))]);
end
end

end

progressbar( h,-1 ); 

strength(k);

perc_sk = prctile(vol_sk,prt);
perc_corr = prctile(vol_corr,prt);



% plot the curves

  figure
  plot(pvals,perc_sk);
    xlabel('Signal size (voxels)'); 
    ylabel('S_o_b_s');
    title([num2str(prt) ' percentile of S_o_b_s values for d.f. = ' num2str(df) ' and ' num2str(N) ' voxels']);

 figure
 plot(pvals,perc_corr);
    xlabel('Signal size (voxels)'); 
    ylabel('r_o_b_s');
    title([num2str(prt) ' percentile of r_o_b_s values for d.f. = ' num2str(df) ' and ' num2str(N) ' voxels']);

fprintf('\n Done. \n\n');