function [file_names, rate_info, rates_info] = rates(vargin)

delete(figure(1)); delete(figure(2)); delete(figure(3));
%cond = input('\n What condition is this? \n [1] uncued tapping \n [2] synch/cont tapping \n [3] rhythm repitition \n [4] plot raw data \n > ');

% Read in all data files

files = spm_select([1 Inf],'.mat$','Select .mat file(s)');
% loop through all files

tap_output = zeros(size(files,1),6);

rates_info = zeros(size(files,1),6);

for k = 1:size(files,1)

file = files(k,:); 

file_name = dir(file);
file_name = struct2cell(file_name);
file_names(k) = file_name(1);

if k == 1
% change to that directory
file = char(file);
namel = size(file,2);

filel = dir(file);
filel = numel(filel.('name'));

new_dir = file(1:(namel - filel));
cd(new_dir);

end

[rate_info] = check_rate(file);

rates_info(k,1) = rate_info.int_mean;
rates_info(k,2) = rate_info.int_mode;
rates_info(k,3) = rate_info.int_min;
rates_info(k,4) = rate_info.int_max;
rates_info(k,5) = rate_info.rate_mean;
rates_info(k,6) = rate_info.rate_mode;
end

file_names = file_names';



