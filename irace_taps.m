function [output] = irace_taps(filenames)

% code to extract timepoints from iRACE outputs
% * should be done on a subject-by-subject basis so that the affected side
% can be appropriately input
%
% rje | 2014.03.27


% select the "touch" files as ground truth
%[filenames, pathname] = uigetfile( {'*touch.txt','TXT-files (*.txt)'},'Select files','MultiSelect', 'on');

% how many files?
% if ischar(filenames)
%    numf = 1;
%    dt = 1;
% elseif iscell(filenames)
%    numf = size(filenames,2);
%    dt = 2;
% end

numf = size(filenames,1);

%% variables to save
fnames = cell(numf,1);
ies_mn = nan(1,numf);
ies_cv = nan(1,numf);
phi_mn = nan(1,numf);
phi_cv = nan(1,numf);
cols   = nan(numf,4); 


%% loop
for f = 1:numf
%    if dt == 1
%        fname = filenames(f,:);
%    elseif dt == 2
%        fname = char(filenames(f));
%    end
%   
%file_touch = [pathname fname]; 

    file_touch = filenames{f};
    [pathname filename] = fileparts(file_touch);   
    fnames{f} = filename;
   
    X_t = load(file_touch);
   
%    % automatically find the corresponding "non-touch" files
%    file_reg = [pathname fname(1:(end-10))];
%    X_r = load(file_reg);
%    
% 
%    % process the regular file
%    % sort to get all touchdown (TD) times in order touchdown = column 1 "0"
%    X_r = sortrows(X_r,[1 3]);
% 
%    % just keep the TD times
%    td_r = find(X_r(:,1) == 0);
%    Xtd_r = X_r(td_r,:);
% 
%    % get the tap series
% 
%    % ------------------------------
%    % 1. left and right separate for phase
%    % column 2: 0 = left, 1 = right
% 
%    indl = find(Xtd_r(:,2) == 0);
%    indr = find(Xtd_r(:,2) == 1);
% 
% 
%     tapsl_r = Xtd_r(indl,3); tapsl_r = tapsl_r(:);
%     tapsr_r = Xtd_r(indr,3); tapsr_r = tapsr_r(:);
     

%% process the touch file
% there will be 5 columns of data: (1) clock; (2) touchdown = 0; (3) ?; (4) x-coord; (5) y-coord
       
   % sort by touchdown, and then by timestamp
   X_t = sortrows(X_t,[2 1]);

   % just keep the TD times
   td_t = find(X_t(:,2) == 0);
   Xtd_t = X_t(td_t,:);

   % determine left vs right; middle = 240 (confirmed)
   xval = Xtd_t(:,4);
   yval = Xtd_t(:,5);

   % get Left vs Right events marked
   button = xval >= 240; % will have L = 0 and R = 1;
   button = button + 1; % now will have L = 1 and R = 2;

    % make the ev_vec format
    ev_vec = [Xtd_t(:,1) button];

    % now we need to exclude a few non-tapping events (continue button etc)
    ev_vec = ev_vec(yval>80,:); % continue button is at the top of the screen, and has y coordinates > 80

    % for L and R tapping, the continue button is actually in a "legal" place
    
    % count number of L and R
    numL = sum(button == 1);
    numR = sum(button == 2);
    
    if numL < 5
        % only count the R button
        ev_vec = ev_vec(button==2,:);
    elseif numR < 5
        % only count the L button
        ev_vec = ev_vec(button==1,:);
    end

    %% outcome measures
    % for now, just skip the rapid files
    if findstr(file_touch,'Rapid') > 0
        % just move on to next file
    else
        outcomes = irace_outcomes(ev_vec);

        ies_mn(f) = outcomes.ies_mn;
        ies_cv(f) = outcomes.ies_cv;
        phi_mn(f) = outcomes.phi_mn;
        phi_cv(f) = outcomes.phi_cv;   
        cols(f,:) = [outcomes.ies_mn outcomes.ies_cv outcomes.phi_mn outcomes.phi_cv];
    end
end % file loop

%% outputs
output.num_files = numf;
output.fnames = fnames;
output.ies_mn = ies_mn;
output.ies_cv = ies_cv;
output.phi_mn = phi_mn;
output.phi_cv = phi_cv;
output.cols = cols;




