function [output pdiff] = swing_time_calcs

% simple function to read in the matrix of swing time data, formatted per
% Hang's code
%
% "pdiff" is the percent differene between L and R midhinges: 100 * (L - R) / (0.5 * (L + R))
%

[filenames, pathname] = uigetfile( {'*.mat','Swing time data'},'Select desired file','MultiSelect', 'off');
cd(pathname)

file = [pathname filenames];
data = load(file);
data = data.Swing; % the full data file can be larger than this, but we just pull out the swing data

pdiff =   cell(size(data,1),1);

for i = 1:size(data,1)
    L = data{i,2}; L = L(:);
    R = data{i,3}; R = R(:);
    
    % percent difference between them
    if isempty(L) || isempty(R)
        % problem
        pdiff{i} = NaN;
    else
        miny = min(min(L),min(R));
        maxy = max(max(L),max(R));
    
        midhingeL = 0.5 * (prctile_nist(L,25) + prctile_nist(L,75));
        midhingeR = 0.5 * (prctile_nist(R,25) + prctile_nist(R,75));
    
        pdiff{i} = 100* (midhingeL - midhingeR) / (0.5 * (midhingeL + midhingeR)); % positive means R > L

        figure(50)
        plot(repmat(1,numel(L),1),L,'.')
        hold on
        plot(1,midhingeL,'+r')
        plot(repmat(2,numel(R),1),R,'.')
        plot(2,midhingeR,'+r')
        hold off
        axis([0 3 miny maxy])
        %pause(0.25)
    end
 
end

output = [data pdiff];
