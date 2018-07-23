function [stats, rvals] = mc_scatter
% mc_scatter
%
% user selects a set of volumes (e.g., T-maps) and perform a desired number
% of 2-combinations of those volumes (e.g., 1000 Monte Carlo iterations), recording the r value for each
% iteration and showing a histogram

% calls: getvol_fxn


% select the files
files = spm_select([1 Inf],'image','Select files for 2-combinations:',[],pwd,'.*');
numf = size(files,1);


% select the masking template
mask=spm_select(1,'image','Select inclusive mask:',[],pwd,'.*');

full_iter = combnk([1:numf],2);
full_count = size(full_iter,1);
tar_iter = input([ '\n ' num2str(full_count) ' iterations possible. How many to perform?: ']);  % 500 seems to be good

if tar_iter > full_count
    fprintf(' Not enough files to do this many 2-combinations! \n');
    return
end

% get the combnk set up

if tar_iter < full_count
    
prm = randperm(full_count); % a random permutation of the total possible pairs
pairs = zeros(tar_iter,2);
for i = 1:tar_iter
    pairs(i,:) = full_iter(prm(i),:);
end

elseif tar_iter == full_count
    % no need to do the permutation if using all of them
pairs = full_iter(1:tar_iter,:);  % just take the first "tar_iter" number of pairs
end

% read in the mask

 vm = spm_read_vols(spm_vol(mask));
 [x1 y1 z1] = size(vm);
 vm2 = vm(:);
 
% create the matrix of all file values for the combnk

fvals_all = zeros(numel(vm2(vm2~=0)),tar_iter);

% read in the file values (files should all be the same size)

for i = 1:numf
   file = files(i,:); 
   [fvals] = getvol_fxn(file,vm2); 
   
   % put the fvals for each volume into the main matrix
   fvals_all(:,i) = fvals;
end

% now we can take all the corr pairs

rvals = zeros(tar_iter,1);

for i = 1:tar_iter
    
   f1 = fvals_all(:,pairs(i,1));
   f2 = fvals_all(:,pairs(i,2));
   
   rvals(i) = corr(f1,f2);
end

figure
hist(rvals,20)

% get the stats for this run

stats.mean = mean(rvals);
stats.std = std(rvals);
stats.prc = prctile(rvals,[95 99 99.5 99.9]);

fprintf(' Finished. (Percentile values are 95, 99, 99.5, 99.9.) \n\n');