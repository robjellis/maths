% Simulating regions of signal in trnd distributions
%
% Just does vector simulations; therefore, assume that FWHM = 0
% sig = proportion of total volume in which signal is put (.0005, .001,
% .0015, .002)

N = input(' Number of voxels: ');
df = input(' Enter the d.f.: ');
iter = input(' Enter the number of iterations: ');
Sobs = input (' Enter the Sobs value: ');

sig = [.0005 .001 .0015 .002];

mag = [4 5 6];
sd = [1 1 1 1];

sk = zeros(iter,length(mag));

vol_temp = zeros(N,2);  

for k = 1:length(sig);
fprintf([' ====== signal size: ' num2str(sig(k)) '\n']);

for j = 1:length(mag);
for i = 1:iter
    
    vol = trnd(df,N,1);
    vol1 = vol;
    vol(1:(ceil(sig(k) *N))) = randn(ceil((sig(k) * N)),sd(j)) + mag(j); % random numbers with sd(j) and mag(j)
    
    sk1(i) = skewness(vol1)
    sk2(i) = skewness(vol)
    
    sk(i,j) = skewness(vol);
    
if mod(i,2)~=0   % if i is odd ...    
vol_temp(1:N,1) = vol';

elseif mod(i,2) == 0 % if i is even ...
vol_temp(1:N,2) = vol';    
end

plot(sk1,sk2)
bef_after(j) = corr(sk1',sk2'); % this is the correlation between all the simulated distributions before and after mag(j) is inserted

end



bef_after

% tally the percentage of iter runs where Sobs > S95

   prt = sum(sk>Sobs) / iter

   
% just find the 5 percentile value

prt5 = prctile(sk,5)

end
end