name=spm_get();
v1=spm_vol(name);
v2=v1.('mat');
origin = [abs(v2(1,4)/v2(1,1)), abs(v2(2,4)/v2(2,2)), abs(v2(3,4)/v2(3,3))]



