function [outvals outvals2] = nullit(m1,m2)

%
% function [outvals] = nullit(m1,m2)
%
% m1 = an A x B matrix
% m2 = an A x B matrix
%
% version = 2011.12.20

vals1 = numel(m1);
vals2 = numel(m2);

if sum(size(m1) - size(m2)) == 0
    % OK
else
    fprintf('\n Warning: vector A and vector B do not have the same length. Terminating.\n\n');
end

% loop ...

outvals = ((abs(m1)<abs(m2)) .* m1) + ((abs(m2)<abs(m1)) .* m2);

% now, zero values where the sign of m1 and m2 does not match

signchk = sign(m1) == sign(m2);

outvals = outvals .* signchk;

% let's see how the distribution of non-zero values changes
outvals2 = outvals(outvals ~= 0);

%figure(70)
%hist(m1,30,'FaceColor','flat')
%hold on
%hist(outvals,30,'FaceColor',[1 0 0])
%hold off



