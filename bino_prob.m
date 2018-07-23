function result = bino_prob(n,p)
%
% n = number of subjects
% p = probability of sucess 
% Function will return the probability of getting >= k sucesses
%

result = []; %clears result variable

%n = input('Total number of subjects:',2,'i');
%p = input('Prob. of success (0-1):',2,'r');

for kk = 1:n
    pk = (1 - (binocdf(kk-1,n,p)));
    result(kk,:) = [kk, pk];
end
fprintf('Result below indicates probability <p> of observing <k> or greater sucesses. \n \n');
fprintf('\t k \t p \n');





