function [output partial_r_all semipartial_r_all] = multreg(inputs,data,regr)

% function [output] = multreg(inputs,data,regr)
%
% inputs: 0 for full simulation, 1 for manual and 2 for SPM
%
% for full simulation:
% data = number of subjects
% regr = number of regressors
%
% for manual:
% data is an [N x 1] matrix of values; N is subjects 
% regr is an [N x M] matrix of values; N is subjects and M is regressors
%
% for SPM:
% (not coded yet)
%
% version = 2012.09.15 (rje)

useit = 0; % rje dummy variable to cut out sections of code that are not needed

if inputs == 0
    % a fully random simulation
    regr = randn(data,regr); % random values
    data = randn(data,1); % random values
    
    
    nsub = size(data,1);
    nreg = size(regr,2); % number of regressors

    
elseif inputs == 1
   % we use "data" and "regr" as existing variables
    nsub = size(data,1);
    nreg = size(regr,2); % number of regressors

    
    % for now, we only do a single vector of data for simplicity
    if size(data,2) ~= 1
       fprintf('\n Note: the program only handles a single vector of [N x 1] subjects data.\n\n');
       return
    end
    
   if nsub == size(regr,1);
       % OK
   else
        fprintf('\n The DATA and REGR matrices do not have the same number of row. Respecify.\n\n');
        return
   end
   
   % Finally, we need to eliminate any cases that have NaN
   xx = 1:nsub; % index
   
   xx = xx(isnan(data)==0); % retain these rows for data and regr

   data = data(xx);
   regr = regr(xx,:);
   nsub = size(data,1); % redefine this now
   
elseif inputs == 2
    
   % this will be for SPM data
   files = spm_select([1 Inf],'image','Select all scans:',[],pwd,'.*');
   mask = spm_select(1,'image','Select the implicit mask:',[],pwd,'.*');
   
    vm = spm_read_vols(spm_vol(mask));
    maskdim = size(vm);
    vm(isnan(vm)) = 0;  % replaces NaNs with zero
    vm = vm ~= 0;       % binarize
    vm = double(vm);
   
   nsub = size(files,1); % number of subjects
   nreg = size(regr,2); % number of regressors

   
   if nsub == size(regr,1);
   % OK
   else
        fprintf('\n The DATA and REGR matrices do not have the same number of row. Respecify.\n\n');
        return
   end

end 

    
% check size



% 3. Options
% note: we either do a true jackknife or we do pure Monte Carlo (random on
% each trial); there is no true bootstrap in this program

opt = input([ ' \n\n ' num2str(nsub) ' subjects retained. Options: \n    [1] Standard analysis \n    [2] Jackknife resampling \n    [3] Monte Carlo permutations / random vectors \n    [4] Subtract vector permutations: ']);

if opt == 1
   niter = 1;
   d = 0;
elseif opt == 2
   
   dd = input([' Choose one: \n   [1] Delete d out of ' num2str(nsub) ' subjects ; \n   [2] Retain j out of ' num2str(nsub) ' subjects: ']);
   if dd == 1
      d = input(' Enter the number of subjects to delete (d): ');
   elseif dd == 2
      d = input(' Enter the number of subjects to retain (j): ');
      d = nsub - d;
   end
   
   niter1 = nchoosek(nsub,d);
   if niter1 > 250000
      niter2 = input(['\n ' num2str(niter1) ' unique subsamples. \n How many *random* subsamples to test?: ' ]);
   else 
      niter2 = input(['\n ' num2str(niter1) ' unique subsamples. \n How many of these to test?: ']);
   end
   
   if niter2 <= niter1
       % OK
       niter = niter2;
   else
       fprintf('\n Cannot evaluate more subsamples than there are possible subsamples. Specify a lower number.\n ')
       return
   end
   
end

if opt == 3
   d = 0;
   datatype = input(' Dependent variable: \n    [1] permutations of input data or \n    [2] randn values: ');
   niter = input('\n Number of Monte Carlo iterations to perform: ');
elseif opt == 4
   d = 0;
   datatype = input(' Permutations of subtracted: \n    [1] input vector or \n    [2] randn values: ');
   niter = input('\n Number of Monte Carlo iterations to perform: ');
   
   if datatype == 1
       
        sub_vec = input([' Enter a [' num2str(nsub) ' x 1] vector, in [ ]: ']); 
   end
   
   if size(sub_vec,1) == size(data,1)
       % OK
   else
       fprintf('\n Warning: the subtraction vector does not have the same size as the data vector. Respecify.\n');
       return
   end
end

% ---------------------
% for the brain data ...

if inputs == 2
    allsub = zeros(maskdim(1), maskdim(2), maskdim(3), nsub);    
    % read in the data
        for i = 1:nsub
            fvol = spm_read_vols(spm_vol(files(i,:)));  
            if sum(size(fvol) - size(mvol)) == 0
               % OK
            else
               % problem
               fprintf('\n Not all files have the same dimensions. ');
            end
            allsub(:,:,:,i) = fvol;   
        end
    
end

fprintf('\n\n Working ...\n\n');

% =====================
% set up variables

% ***** for regression module to work in matlab, we must have a vector of 1s in
% the regressor matrix; i.e., the intercept
% see: http://www.mathworks.com/help/toolbox/stats/regress.html

regrmod = regr;
regrmod(1:nsub,nreg+1) = ones(nsub,1); 

if inputs == 0 || inputs == 1
    zr = zeros(niter,nreg); % between each regressor and the data (not between different regressors)
    pr = zeros(niter,nreg); 
    sr = zeros(niter,nreg); % individual runs as appropriate

    R2 = zeros(niter,1); % the actual R2 for the model iterations
    adjR2 = zeros(niter,1); % adjusted R2 (Ezekiel formula)
    Fval = zeros(niter,1); % the associated F-statistic


elseif inputs == 2
    % just get the R2 volume at this point
    R2vol = zeros(maskdim);
end

if inputs == 0 || inputs == 1
    
    if opt == 2 || opt == 3 || opt == 4
        progressbar
    end
    
    if opt == 2
        anal = ['delete-' num2str(d) ' jackknife'];

        if niter1 <= 250000 % ALL possible sub-samples
 
            xin = nchoosek(1:nsub,nsub-d); % all the possible combinations of nsub-d subjects

            % we may only want a subsample of this
            rx = randperm(niter1);
            rx = rx(1:niter); % just the first niter of the total possible, after randomization

            xin = xin(rx,:);
        else % there are too many possible subsamples, so we define a random set here
            xin = zeros(niter2,nsub-d);
            for xx = 1:niter2
                xtemp = randperm(nsub); % random permutation of all subjects
                xtemp = xtemp(1:nsub-d); % just the cases we need
                xin(xx,:) = xtemp;
            end
        end
    end
    
    % ========================
    % the loop
    for i = 1:niter
        if opt == 1 % original data, no mods necessary
           
           data2 = data;
           regrmod2 = regrmod;
           zc = corr(data2,regrmod2);
           zc = zc(1:end-1);
              % zero-order correlations
              zr(i,:) = zc;
        elseif opt == 2 % bootstrap

            %pause(0.01) % Do something important
            progressbar(i/niter) % Update progress bar

            
            xtrial = xin(i,:); % the cases to keep for this iteration
            data2 = data(xtrial); % the subset of the data cases
            regrmod2 = regrmod(xtrial,:); % the subset of the regressor values
            
            zc = corr(data2,regrmod2);
            zc = zc(1:end-1);
            zr(i,:) = zc;

        elseif opt == 3 || opt == 4

            %pause(0.01) % Do something important
            progressbar(i/niter) % Update progress bar

            %if d > 0
                % get the subset
                %x = randperm(nsub);
                %x = x(1:nsub-d); % a UNIQUE *random* subset on every iteration
            %elseif d == 0
                %x = 1:nsub; % the full set; permutation happens next
            %end
            
            % note: below, data and regr are still "matched" in the correct
            % order (i.e., a jackknife)
            %data3 = data(x); % retain only those subjects in x
            %regrmod2 = regrmod(x,:);

           % for MC and subtraction permutations (not delete-d jackknife)
           if opt == 3
                if datatype == 1 % for MC
                    anal = 'M.C. using permutations';
                    xperm = randperm(nsub); % a UNIQUE random permutation of integers 1:nsub for each iteration
                    data2 = data(xperm);    % a reordering of the data
                    regrmod2 = regrmod(xperm,:);
                    zc = corr(data2,regrmod2);
                    zc = zc(1:end-1);
                    zr(i,:) = zc;
                    
                elseif datatype == 2 % for MC
                    anal = 'M.C. using randn values';
                    data2 = randn(nsub,1);
                    
                    % new regressors on every iteration as well
                    regrmod2 = randn(nsub,nreg);
                    
                    % zero-order correlations
                    zr(i,:) = corr(data2,regrmod2);
                    
                    % intercept column
                    regrmod2(1:nsub,nreg+1) = ones(nsub,1); 
                end
                
           elseif opt == 4
                
                if datatype == 1 % for the subtraction
                    anal = 'Data minus vector subtractions';
                    xperm = randperm(nsub); % a UNIQUE random permutation of integers 1:nsub-d for each iteration
                    data2 = data - sub_vec(xperm);
                elseif datatype == 2
                    anal = 'Data minus randn subtractions';
                    data2 = randn(nsub,1);
                    data2 = data - data2;
                end
           end
            
        end
        
        % ----------------------
        % run the full regression
        [b,bint,resid,rint,stats] = regress(data2,regrmod2);
        
        Fval(i) = stats(2); % the associated F-value for this model

        % the R2 is the square of the correlation between 
        % (1) the original data (y-values) and
        % (2) those values minus the residuals from the full model
        % (residualized x-values)
        %
        % it is equal to stats(1) returned above
        
        if niter == 1
            %figure(20)
            %plot(data2,data2-resid,'.')
            %xlabel('data'); ylabel('data - residuals');
        end
        
        R2(i) = (corr(data2,data2-resid))^2;
        adjR2(i) = 1 - (1-R2(i))*(nsub-1)/(nsub - nreg - 1); % the Ezekiel (1930) formula widely used by SPSS, Statistica, etc.
                                                             % Other formulas are available, as discussed in Leach (2007)
        % rje: if adjR2 is < 0, just make it zero; this is what Leach does
        % thus negative values won't hurt the mean of the distribution
        if adjR2(i) < 0
           adjR2(i) = 0;
        end
        
      % ----------------------  
      % semipartial correlation analysis
        for m = 1:nreg % for each of the regressors in succession ...
           y = 1:nreg+1; % need to include the constant vector
           y = y(y~=m);  % the other regressors
           
           % now we get the semipartial r by correlating 
           % (1) the original data (y-values) and
           % (2) residualized x-values (adjusted for IVs in regression model)
           
           regrmod3 = regrmod2(:,y); % still includes the final constant vector
           [b,bint,resid] = regress(regrmod2(:,m),regrmod3);
          
           sr(i,m) = corr(data2,resid);

        end
      
      % ----------------------
      % partial correlation analysis 
       for m = 1:nreg
       y = 1:nreg;
       y = y(y~=m);

       X = [data2 regrmod2(:,m)];
       Y = regrmod2(:,y); % don't include the regressor just put into X

       pc = partialcorr(X,Y);
       pr(i,m) = pc(1,2); % just the correlation value

       end
       

    end % iteration loop
       
       % how many partial r (and semipartial) r values are significnat at .05?
       rpvals = r2p(pr,'b',nsub-d,nreg); % these are now p-values       
       sig_rpvals = sum(rpvals < .05);
       
       % now turn this into a percent
       numr = size(pr,1);
       partial_r_persig = 100 * (sig_rpvals / numr);
       
       
       % what is the correlation between model R and (1) max partial r and 
       % and (2) mean semipartial r?
       
       %model_R = sqrt(R2);
       partial_r2 = power(pr,2); % these are the r2 values
       partial_r_all = pr(:);
       semipartial_r2 = power(sr,2); % these are the r2 values
       semipartial_r_all = sr(:);
       zr_all = zr(:);
       
       % get the critical r-values; default is p = .05 with 2 tails
       rc = rcrit(nsub-d,nreg,.05,2);
       zc = r2z(rc); % fisher z-values
       
       % relationship between zero-order, partial, and semipartial
       if inputs == 0
           figure(21)  
           subplot(1,2,1)    
           plot(zr_all,partial_r_all,'.')
           xlabel('zero-order r'); ylabel('partial r');
           axis([-1 1 -1 1])      

           subplot(1,2,2)    
           plot(zr_all,semipartial_r_all,'.')
           xlabel('semipartial r'); ylabel('partial r');
           axis([-1 1 -1 1]) 

           %plot(max(pr,[],2),model_R,'.')
           %xlabel('max partial r'); ylabel('model R');


           % what is the relationship between partial r and semipartial r?
           figure(22)
           plot(partial_r_all,semipartial_r_all,'.')
           xlabel('partial r'); ylabel('semipartial r');
           axis([-1 1 -1 1])

           % what are the distributions of these values?
           %xmin = min(min(partial_r_all),min(semipartial_r_all));
           %xmax = max(max(partial_r_all),max(semipartial_r_all));
           xvals = linspace(-1,1,50); 

           y1 = histc(partial_r_all,xvals);
           y2 = histc(semipartial_r_all,xvals);
           y3 = histc(zr_all,xvals);

           figure(23)
           plot(xvals,y3,'g','LineWidth',3,'DisplayName','zero-order r')
           hold on
           plot(xvals,y1,'b','LineWidth',2,'DisplayName','partial r')
           plot(xvals,y2,'r','LineWidth',1,'DisplayName','semipartial r')
           xlabel('r-value'); ylabel('count');
           legend('show', 'Location','NorthWest');
           hold off
       
              
       % convert to Fisher z-values
       partial_z = 0.5 .* log((1+pr)./(1-pr));
       semipartial_z = 0.5 .* log((1+sr)./(1-sr));
       
       
       xrange = [0 nreg+1];
       rpos = [rc rc];
       rneg = -1*rpos;
       zpos = [zc zc];
       zneg = -1*zpos;
       
       % zero-order r, partial r, and semipartial r
       figure(24)
       subplot(1,3,1)
       plot(xrange,rpos,'g','LineWidth',2)
       hold on
       plot(xrange,rneg,'g','LineWidth',2)
       boxplot(zr,'notch','off','whisker',.953,'symbol','r')
       title('Zero-order r')
       
       subplot(1,3,2)
       plot(xrange,rpos,'g','LineWidth',2)
       hold on
       plot(xrange,rneg,'g','LineWidth',2)
       boxplot(pr,'notch','off','whisker',.953,'symbol','r')
       title('Partial r')
       
       subplot(1,3,3)
       plot(xrange,rpos,'g','LineWidth',2)
       hold on
       plot(xrange,rneg,'g','LineWidth',2)      
       boxplot(sr,'notch','off','whisker',.953,'symbol','r')
       title('Semipartial r')
       
       end % for plotting
       
       % rje calculated .953 to give whiskers at .025 and .975, yielding a
       % 95% interval. Default is to ignore "outliers"; add a + after 'r'
       % to put them back in
       
       % partial z and semipartial z
       
       if useit == 1
       figure(25)
       subplot(1,2,1)
       plot(xrange,zpos,'g','LineWidth',2)
       hold on
       plot(xrange,zneg,'g','LineWidth',2)
       boxplot(partial_z,'notch','off','whisker',.953,'symbol','r')
       title('Partial r2z')
       
       subplot(1,2,2)
       plot(xrange,zpos,'g','LineWidth',2)
       hold on
       plot(xrange,zneg,'g','LineWidth',2)      
       boxplot(semipartial_z,'notch','off','whisker',.953,'symbol','r')       
       title('Semipartial r2z')
       end
       

       
       % what is the relationship among partial r2 values across the set of
       % iterations?
       partial_r2_self_corr = corr(partial_r2);
       
       % hermite spline of the CDF so that we can have more exact values
       %xher = 0:.0001:1;
       %r_cdf_her = pchip(xvals,y1,xher); 
       %sr_cdf_her = pchip(xvals,y2,xher);
       

    
elseif inputs == 2
    
   % run the full regression
   
   
   % apply the mask
   
   vol = vol .* vm;
   
 % write the R2 volume
 v1n.('fname') = strcat(file1, '_Bin_',num2str(vec1),'.nii');    % will write to the "file1" directory
 v1n.('pinfo') = [0; 0; 0];                    % reset this field so that scale factors are not preserved
 v1n = spm_write_vol(v1n,vol);
   
end
  


perc = [2.5 97.5]; % to create a 95% CI

% summary statistics
% note: F-distributions are one-tailed, since F-values cannot be negative




if niter == 1
    
    Fcrit_P05 = finv(1-.05,nreg,(nsub-nreg-1));
    Fcrit_P01 = finv(1-.01,nreg,(nsub-nreg-1));
    Fcrit_P005 = finv(1-.005,nreg,(nsub-nreg-1));
    Fcrit_P001 = finv(1-.001,nreg,(nsub-nreg-1));

   output.R2 = R2;
   output.adjR2 = adjR2;
   
   output.Fval = Fval;
   output.Fcrit_P05 = Fcrit_P05;
   output.Fcrit_P01 = Fcrit_P01;
   output.Fcrit_P005 = Fcrit_P005;
   output.Fcrit_P001 = Fcrit_P001;
   
   output.zr = zr;
   output.pr = pr;
   output.sr = sr;
   output.P_vals = rpvals;
   
elseif niter > 1
   output.analysis = anal;
   output.num_cases = nsub - d;
   output.iterations = niter;
   output.percentile = perc;
   output.r_crit_2_tail = rc;
   output.x0 = '               ';
   output.R2_med = median(R2);
   output.R2_sd = std(R2);
   output.R2_perc = prctile(R2,perc);
   output.x1 = '               ';
   output.adjR2_med = median(adjR2);
   output.adjR2_sd = std(adjR2);
   output.adjR2_perc = prctile(adjR2,perc);
   %output.x2 = '               ';
   %output.Fval_mn = mean(Fval); 
   %output.Fval_sd = std(Fval);
   %output.Fval_perc = prctile(Fval,perc);
   %output.Fcrit_P05 = Fcrit_P05;
   %output.Fcrit_P01 = Fcrit_P01;
   %output.Fcrit_P005 = Fcrit_P005;
   %output.Fcrit_P001 = Fcrit_P001;

   output.x2 = '               ';
   output.zr_med = median(zr);
   output.zr_sd = std(zr);
   output.zr_perc = prctile(zr,perc);
   output.zr_all = zr;

   output.x3 = '               ';
   output.partial_r_med = median(pr);
   output.partial_r_sd = std(pr);
   %output.partial_r_p_mn = mean(rpvals);
   %output.partial_r_perc_hi = prctile(pr,perc(2));
   output.partial_r_perc = prctile(pr,perc);
   output.partial_r_persig = partial_r_persig;
   output.partial_r_all = pr;
   %output.partial_z_mn = mean(partial_z);
   %output.partial_z_sd = std(partial_z);
   %output.partial_z_perc = prctile(partial_z,perc);
   %output.partial_r_self_corr = partial_r_self_corr;
   output.x4 = '               ';
   output.semipartial_r_med = median(sr);
   output.semipartial_r_sd = std(sr);
   %output.semipartial_r_perc_hi = prctile(sr,perc(2));
   output.semipartial_r_perc = prctile(sr,perc);
   %output.semipartial_z_mn = mean(semipartial_z);
   %output.semipartial_z_sd = std(semipartial_z); 
   %output.semipartial_z_perc = prctile(semipartial_z,perc);
   output.semipartial_r_all = sr;
   output.corr_R2__partial_r2 = corr(R2,max(partial_r2,[],2));
   output.corr_R2__semipartial_r2 = corr(R2,max(semipartial_r2,[],2));
   output.corr_partial_r2__semipartial_r2 = corr(partial_r2(:),semipartial_r2(:));
       
     
end




