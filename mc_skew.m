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

for k = 1:length(strength)

    h = progressbar( h,(1/length(strength)));                 
    
    if ~gui_active
        break;
    end  
    
for j = 1:length(sig)
    
for i = 1:iter
    
  vol1 = trnd(df,N,1);
  vol1(1:sig(j)) = randn(sig(j),1) + strength(k);
  
  vol2 = trnd(df,N,1);
  vol2(1:sig(j)) = randn(sig(j),1) + strength(k);
  
  vol_mn(i, j) = mean(vol1);
  vol_sk(i, j) = skewness(vol1);
  vol_corr(i, j) = corr(vol1,vol2);

if makecorr == 'y' && i == iter && j == length(sig)
   figure
   plot(vol1,vol2,'.');
   xlabel('Values from vol1'); 
   ylabel('Values from vol2');
   title([num2str(sig(j)) ' voxels with M = ' num2str(strength(k))]);
end
end

end
strength(k);

if k == 1
a_skew_vol = prctile(vol_sk,prt);
a_corr_vol = prctile(vol_corr,prt);
clear vol_mn vol_sk vol_corr
elseif k == 2
b_skew_vol = prctile(vol_sk,prt);
b_corr_vol = prctile(vol_corr,prt);
clear vol_mn vol_sk vol_corr
elseif k == 3
c_skew_vol = prctile(vol_sk,prt);
c_corr_vol = prctile(vol_corr,prt);    
clear vol_mn vol_sk vol_corr
elseif k == 4
d_skew_vol = prctile(vol_sk,prt);
d_corr_vol = prctile(vol_corr,prt);    
clear vol_mn vol_sk vol_corr
elseif k == 5
e_skew_vol = prctile(vol_sk,prt);
e_corr_vol = prctile(vol_corr,prt);    
clear vol_mn vol_sk vol_corr
elseif k == 6
f_skew_vol = prctile(vol_sk,prt);
f_corr_vol = prctile(vol_corr,prt);

end
end

progressbar( h,-1 ); 

% plot the skewness curves

  figure
  plot1 = plot(sig,a_skew_vol);
  hold on
  plot2 = plot(sig,b_skew_vol);
  plot3 = plot(sig,c_skew_vol);
  plot4 = plot(sig,d_skew_vol);
  plot5 = plot(sig,e_skew_vol);
  plot6 = plot(sig,f_skew_vol);
    xlabel('Signal size (voxels)'); 
    ylabel('S_o_b_s');
    title([num2str(prt) ' percentile of S_o_b_s values for d.f. = ' num2str(df) ' and ' num2str(N) ' voxels']);
    
    legend1 = legend('show');
    set(legend1,'Location','NorthWest');

    set(plot1,'DisplayName',num2str(strength(1)),'Color',[1 0 0],'Marker','.');
    set(plot2,'DisplayName',num2str(strength(2)),'Color',[0.8706 0.4902 0],'Marker','.');
    set(plot3,'DisplayName',num2str(strength(3)),'Color',[0 0.498 0],'Marker','.');
    set(plot4,'DisplayName',num2str(strength(4)),'Color',[0 0.749 0.749],'Marker','.');
    set(plot5,'DisplayName',num2str(strength(5)),'Color',[0 0 1],'Marker','.');    
    set(plot6,'DisplayName',num2str(strength(6)),'Color',[0.4784 0.06275 0.8941],'Marker','.');

  hold off
  
  
% plot the r-value curves  
  figure
  plot1 = plot(sig,a_corr_vol);
  hold on
  plot2 = plot(sig,b_corr_vol);
  plot3 = plot(sig,c_corr_vol);
  plot4 = plot(sig,d_corr_vol);
  plot5 = plot(sig,e_corr_vol);
  plot6 = plot(sig,f_corr_vol);
    xlabel('Signal size (voxels)'); 
    ylabel('r_o_b_s');
    title([num2str(prt) ' percentile of r_o_b_s values for d.f. = ' num2str(df) ' and ' num2str(N) ' voxels']);
    
    legend1 = legend('show');
    set(legend1,'Location','NorthWest');

    set(plot1,'DisplayName',num2str(strength(1)),'Color',[1 0 0],'Marker','.');
    set(plot2,'DisplayName',num2str(strength(2)),'Color',[0.8706 0.4902 0],'Marker','.');
    set(plot3,'DisplayName',num2str(strength(3)),'Color',[0 0.498 0],'Marker','.');
    set(plot4,'DisplayName',num2str(strength(4)),'Color',[0 0.749 0.749],'Marker','.');
    set(plot5,'DisplayName',num2str(strength(5)),'Color',[0 0 1],'Marker','.');    
    set(plot6,'DisplayName',num2str(strength(6)),'Color',[0.4784 0.06275 0.8941],'Marker','.');
    
  hold off

fprintf('\n Done. \n\n');