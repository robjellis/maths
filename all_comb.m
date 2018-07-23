function out = all_comb(varargin)

% enumerate all combinations of the data
% for use in the StarPower project
%
% RJE | 10 Feb 2018

% how many?
nitem = numel(varargin);

if nitem == 2
	% write code if needed
elseif nitem == 5
	i1 = varargin{1}; n1 = numel(i1);
	i2 = varargin{2}; n2 = numel(i2);
	i3 = varargin{3}; n3 = numel(i3);
	i4 = varargin{4}; n4 = numel(i4);
	i5 = varargin{5}; n5 = numel(i5);
	
	% how many total?
	nposs = n1*n2*n3*n4*n5;
	
	labels = nan(nposs,nitem); % matrix for later use (ANOVA etc)
	combs  = nan(nposs,nitem);
	
	ctr = 1;
	
	for p = 1:n1
		for q = 1:n2
			for r = 1:n3
				for s = 1:n4
					for t = 1:n5
						
						labels(ctr,:) = [p q r s t];
						combs(ctr,:) = [i1(p) i2(q) i3(r) i4(s) i5(t)];
						ctr = ctr + 1;
					
					end
				end
			end
		end
	end
	
end

out.levels = labels;
out.combs  = combs;