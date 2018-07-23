function [output] = robustreg(inputs,data_all,regr)

% function [output] = robustreg(inputs,data,regr)
%
% inputs: 0 for full simulation, 1 for manual and 2 for SPM
%
% for full simulation:
% data = number of subjects
% regr = number of regressors
%
% for manual:
% data is an [N x K] matrix of values; N is subjects and K is the number of unique dependent variables 
% regr is an [N x M] matrix of values; N is subjects and M is regressors
%
% for SPM:
% (not coded yet)
%
% version = 2012.09.24 (rje)

useit = 0; % rje dummy variable to cut out sections of code that are not needed

if inputs == 0
    % a fully random simulation
    regr = randn(data_all,regr); % random values
    data = randn(data_all,1); % random values
    
    
    nsub = size(data,1);
    nreg = size(regr,2); % number of regressors

    
elseif inputs == 1
   % we use "data" and "regr" as existing variables
    nsub = size(data_all,1);
    nvar = size(data_all,2); % there will be a loop for multiple data columns
    nreg = size(regr,2); % number of regressors
    
    % make a copy of the original regressor matrix
    regr_orig = regr;
 
   if nsub == size(regr,1);
       % OK
   else
        fprintf('\n The DATA and REGR matrices do not have the same number of row. Respecify.\n\n');
        return
   end
   

   
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
        fprintf('\n The DATA and REGR matrices do not have the same number of rows. Respecify.\n\n');
        return
   end

end 

mc = input(' Mean center each IV?: [1] yes; [2] no: ');    

% 3. Options
% note: we either do a true jackknife or we do pure Monte Carlo (random on
% each trial); there is no true bootstrap in this program

opt = input([ ' Robust regression options: \n    [1] Standard (one-pass) analysis \n    [2] Jackknife resampling \n    [3] Monte Carlo permutations / random vectors \n    [4] Subtract vector permutations: ']);

if opt == 1
   anal = 'Standard robust (all valid cases)';
   niter = 1;
   d = 0;
elseif opt == 2
   if nvar == 1
       % how many valid cases do we have?
       ncase = sum(isnan(data_all)==0);
       fprintf(['\n There are ' num2str(ncase) ' valid cases.\n']);  
   end
   
      % what is the minimum number of valid subjects across the variables?
   all_valid = sum(isnan(data_all)==0)
   min_valid = min(all_valid);
   
   dd = input([' Choose one: \n   [1] Delete d out of N subjects; \n   [2] Retain j out of N subjects: ']);
   if dd == 1
      d = input(' Enter the number of subjects to delete (d): ');
   elseif dd == 2
      d = input(' Enter the number of subjects to retain (j): ');
      d = min_valid - d;
   end
   
   niter1 = nchoosek(min_valid,d);
   if niter1 > 250000
      niter2 = input(['\n ' num2str(niter1) ' unique subsamples across all DVs. \n How many *random* subsamples to test?: ' ]);
   else 
      niter2 = input(['\n ' num2str(niter1) ' unique subsamples across all DVs. \n How many of these to test?: ']);
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
   
   if size(sub_vec,1) == size(data_all,1)
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

% note: robustfit does not require that the intercept vector be added
% see: http://www.mathworks.com/help/stats/robustfit.html

% regular multiple regression in matlab DOES require this
% see: http://www.mathworks.com/help/toolbox/stats/regress.html

% set up the main output DV variables

   nsubs = zeros(nvar,1);
   
   if niter == 1
       iB_all = zeros(nvar,1);
       iSE_all = zeros(nvar,1);
       it_all = zeros(nvar,1);    % regression intercept t-values
       ip_all = zeros(nvar,1);    % associated p-value
       
       pB_all = zeros(nvar,nreg);
       pSE_all = zeros(nvar,nreg);
       pr_all = zeros(nvar,nreg);
       pt_all = zeros(nvar,nreg);
       pz_all = zeros(nvar,nreg);
       pp_all = zeros(nvar,nreg);
   elseif niter > 1
       pr_med = zeros(nvar,nreg);
       pr_sd = zeros(nvar,nreg);
       pr_perc_lo = zeros(nvar,nreg);
       pr_perc_hi = zeros(nvar,nreg);
       
       pp_med = zeros(nvar,nreg);
       pp_sig = zeros(nvar,nreg);

   end


% do the DV loop

for dv = 1:nvar
    
   % need to eliminate any cases that have NaN
   data = data_all(:,dv);
   %size(data)
   
   xx = 1:size(data,1); % index
   xx = xx(:);
   %size(isnan(data))
   yy = xx(isnan(data)==0); % retain these rows for data and regr
   
   data = data(yy);
   regr = regr_orig(yy,:); % take from the ORIGINAL regressor matrix
   nsub = size(data,1); % redefine this now
   nsubs(dv) = nsub;
   
   % mean center the subset regressor matrix?
   
   if mc == 1
        % mean center it (code by rje)
        mndata = mean(regr);
        nn = size(regr,1);
        regr = regr - (ones(nn,1) * mndata);
        
   elseif mc == 2
        % do nothing
   end
   
if inputs == 0 || inputs == 1
    
    pt = zeros(niter,nreg); % regression slope t-values
    pz = zeros(niter,nreg); % conv. to Fisher Z-values
    pr = zeros(niter,nreg); % regression slope partial r-values
    pp = zeros(niter,nreg); % p-values
    %sr = zeros(niter,nreg); % individual runs as appropriate

    %R2 = zeros(niter,1); % the actual R2 for the model iterations
    %adjR2 = zeros(niter,1); % adjusted R2 (Ezekiel formula)
    %Fval = zeros(niter,1); % the associated F-statistic


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
           regr2 = regr;
           %zc = corr(data2,regr2);
           %zc = zc(1:end-1);
              % zero-order correlations
              %zr(i,:) = zc;
        elseif opt == 2 % bootstrap

            %pause(0.01) % Do something important
            progressbar(i/niter) % Update progress bar

            
            xtrial = xin(i,:); % the cases to keep for this iteration
            data2 = data(xtrial); % the subset of the data cases
            regr2 = regr(xtrial,:); % the subset of the regressor values
            
            %zc = corr(data2,regr2);
            %zc = zc(1:end-1);
            %zr(i,:) = zc;

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
            %regr2 = regr(x,:);

           % for MC and subtraction permutations (not delete-d jackknife)
           if opt == 3
                if datatype == 1 % for MC
                    anal = 'M.C. using permutations';
                    xperm = randperm(nsub); % a UNIQUE random permutation of integers 1:nsub for each iteration
                    data2 = data(xperm);    % a reordering of the data
                    regr2 = regr(xperm,:);
                    %zc = corr(data2,regr2);
                    %zc = zc(1:end-1);
                    %zr(i,:) = zc;
                    
                elseif datatype == 2 % for MC
                    anal = 'M.C. using randn values';
                    data2 = randn(nsub,1);
                    
                    % new regressors on every iteration as well
                    regr2 = randn(nsub,nreg);
                     
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
        % run the robust regression
        [b stats] = robustfit(regr2,data2);
        
        se = stats.se; % SE estimates of B values
        ppt = stats.t'; % t-values
        ppp = stats.p'; % p-values
        
        % intercepts
        iB = b(1);
        iSE = se(1);
        it = ppt(1);
        ip = ppp(1);


        % slopes (partial r-values)
        pB = b(2:nreg+1);
        pB = pB';
        pSE = se(2:nreg+1);
        pSE = pSE';
        
           ppt = ppt(2:nreg+1);
        pt(i,:) = ppt;
  
           ppp = ppp(2:nreg+1);
        pp(i,:) = ppp;
        
           ppr = ppt ./ sqrt(nsub - d - nreg - 1 + ppt.^2); % Standard formula for converting r-value to t-value
        pr(i,:) = ppr;
        
        % now convert to Fisher z    
        pz(i,:) = 0.5 * log((1+ppr) ./(1-ppr));

    end % iteration loop
       
       % how many partial r are significant at .05?  
       sig_pvals = sum(pp < .05);
       
       % now turn this into a percent
       percent_sig = 100 * (sig_pvals / niter);
       
       
       % what is the correlation between model R and (1) max partial r and 
       % and (2) mean semipartial r?
       
       %model_R = sqrt(R2);
       %partial_r2 = power(pr,2); % these are the r2 values
       %partial_r_all = pr(:);
       %semipartial_r2 = power(sr,2); % these are the r2 values
       %semipartial_r_all = sr(:);
       %zr_all = zr(:);
       
       % get the critical r-values; default is p = .05 with 2 tails
       rc = rcrit(nsub-d,nreg,.05,2);
       zc = 0.5 * log((1+rc) ./(1-rc));

       
       % rje calculated .953 to give whiskers at .025 and .975, yielding a
       % 95% interval. Default is to ignore "outliers"; add a + after 'r'
       % to put them back in
       
       % partial z and semipartial z
       xrange = [0 nreg+1];
       rpos = [rc rc];
       rneg = -1*rpos;
       zpos = [zc zc];
       zneg = -1*zpos;
       
       if niter > 1
           figure(25)
           subplot(1,2,1)
           plot(xrange,rpos,'g','LineWidth',2)
           hold on
           plot(xrange,rneg,'g','LineWidth',2)
           boxplot(pr,'notch','off','whisker',.953,'symbol','r')
           title('Partial r')

           subplot(1,2,2)
           plot(xrange,zpos,'g','LineWidth',2)
           hold on
           plot(xrange,zneg,'g','LineWidth',2)
           boxplot(pz,'notch','off','whisker',.953,'symbol','r')
           title('Fisher z')
 
       end
       

       
       % what is the relationship among partial r2 values across the set of
       % iterations?
       %partial_r2_self_corr = corr(partial_r2);
       
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

   iB_all(dv,1) = iB;
   iSE_all(dv,1) = iSE;
   it_all(dv,1) = it;
   ip_all(dv,1) = ip;
   
   pB_all(dv,:) = pB;
   pSE_all(dv,:) = pSE;
   pr_all(dv,:) = pr;
   pt_all(dv,:) = pt;
   pz_all(dv,:) = pz;
   pp_all(dv,:) = pp;
   
elseif niter > 1

   %output.rcrit_2tail = rc;
   %output.zcrit_2tail = zc;
   %median(pr)   
   %size(pr_med)
   pr_med(dv,:) = median(pr);
   pr_sd(dv,:) = std(pr);
   pr_perc_lo(dv,:) = prctile(pr,perc(1));
   pr_perc_hi(dv,:) = prctile(pr,perc(2));
   
   % get the p-value for the median r-value
   mm = median(pr);
   tt = mm .* sqrt((nsub-d-nreg-1)./(1- (mm.^2)));
   pp_med(dv,:) = 2*tcdf(-abs(tt),nsub-d-nreg-1); % two-tailed
   
   pp_sig(dv,:) = percent_sig;
   
   %output.pr_all = pr;
   
   %output.pz_med = median(pz);
   %output.pz_sd = std(pz);
   %output.pz_perc_lo = prctile(pz,perc(1));
   %output.pz_perc_hi = prctile(pz,perc(2));
   %output.pz_all = pz;
   
   %output.percent_sig = percent_sig;
  
end % iter choice
   
end % dv loop

output.analysis = anal;
output.num_cases = nsubs - d;
output.df = nsubs - d - nreg - 1;
output.num_reg = nreg;
output.num_iter = niter;


if niter == 1
    output.iB = iB_all;
    output.iSE = iSE_all;
    output.it = it_all;
    output.ip = ip_all;
    
    output.pB = pB_all;
    output.pSE = pSE_all;
    output.pr = pr_all;
    output.pt = pt_all;
    output.pz = pz_all;
    output.pp = pp_all;
    
    
elseif niter > 1
    
   output.pr_med = pr_med;
   output.pr_sd = pr_sd;
   output.prctile = perc;
   output.pr_perc_lo = pr_perc_lo;
   output.pr_perc_hi = pr_perc_hi;
   output.pp_med = pp_med;
   output.pp_sig = pp_sig;
   
   if nvar == 1
      output.pr_all = pr; % for further testing if desired 
   end
end



