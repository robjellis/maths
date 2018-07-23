function [files] = copyfiles

% [files] = copyfiles
%
% copy a list of file paths to a target directory
% works for SPM5 and SPM8
%
% ** note: program assumes all to-be-copied files have a *.xxx extension 
% (i.e., not configured to work with *.nii.gz)
%
% version = 2012.07.07


%% hello

loc = which('copyfiles'); % note: must have a unique name and NOT an internal variable name
                      % in order to have the function called "beats", the
                      % internal variable must have a different name; i.e.,
                      % "beatsvar"
file_info = dir(loc);
save_date = file_info.date;

fprintf(['\n\n || Version: ' save_date '\n || http://robjellis.net \n']);


%%
typ = input('\n Operation: \n [1] Copy N .img/.nii files to a single directory; \n [2] Copy N .img/.nii files to N directories; \n [3] Copy N .mat/.txt files to a single directory \n [4] Make N directories: ');

if typ == 1
    files = spm_select([1 Inf],'image','Select the to-be-copied files or file paths:',[],pwd,'.*');
    numf = size(files,1);
    newdir = spm_select(1,'dir','Select target directory for copied files:');
elseif typ == 2
    files = spm_select([1 Inf],'image','Select the to-be-copied files or file paths:',[],pwd,'.*');
    numf = size(files,1);
    alltar = spm_select([1 numf],'dir','Select the full list of target directories:');
    changename = input(' Rename all copied files to a single file name? [1] yes; [2] no: ');
    % use [1] when setting up a robust regression
    if changename == 1
        newfilename = input(' Enter the new file name to apply to ALL files (no quotes): ','s');
    elseif changename == 2
        % will copy the existing file name
    end
elseif typ == 3
    files = spm_select([1 Inf],'mat','Select the to-be-copied files or file paths:',[],pwd,'.*');
    numf = size(files,1);
    newdir = spm_select(1,'dir','Select target directory for copied files:');
elseif typ == 4
    % prompt for what you need
    tardir = spm_select(1,'dir','Select the parent directory for the new directories:');
    numf = input(' How many new directories to make? (e.g., 50): ');
    dirprefix = input(' Prefix for directories (no quotes needed): ','s');
end


% ---
% set up the counter ...
counter = cell(numf,1);   % will always be _ _ _ (numbers up to 999)

% counter will work either for moving N files or separate files to N
% directories
for q = 1:numf
    
       qq = num2str(q);   % will increment for every unique file   

       if q < 10
          counter{q} = strcat('00',qq);
       elseif q < 100
          counter{q} = strcat('0',qq); 
       elseif q > 99
          counter{q} = strcat(qq); 
       end  

end

fprintf('\n Working ...');

% now we switch 

if typ == 4
   for f = 1:numf 
           qq = counter(f);
           dirnamefin = char(strcat(tardir,dirprefix,'_',char(qq)));
           mkdir(dirnamefin);        
   end
   
elseif typ == 1 || typ == 3
    % may as well switch to the sole target directory
    newdir = char(newdir);
    cd(newdir);

    for f = 1:numf

        file = files(f,:);
        file = strcat(file); % to get rid of trailing spaces
        fend = numel(file);

        % get file name etc - are we using linux or windows?
        if isempty(findstr(file,'/')) == 1
            % then we are using windows
            len = max(strfind(file,'\'));

        elseif isempty(findstr(file,'\')) == 1
            % then we are in linux
            len = max(strfind(file,'/'));
        end

        % do we have a ",1" tag?
           if strcmp(file(fend-1:fend),',1')
              file = file(1:fend-2);
              fend = fend - 2; % respecify this
           end

        % if dealing with .hdr/.img, prepare the .img
        if strcmp(file(fend-3:fend),'.img')
            %fname = file(len+1:fend-3)  % the actual file name after the last divider
            c = 0;
        else % any other format, including .nii, .xls, etc
            c = 1;    
        end   

        fname = file(len+1:fend);  % the actual file name after the last divider

        if c == 0 % .hdr .img format

           file2 = file(1:fend-4);  

           file2img = strcat(file2,'.img');
           file2hdr = strcat(file2,'.hdr');

           fname = file(len+1:fend-4);

           fnameI = strcat(fname,'.img');
           fnameH = strcat(fname,'.hdr');

           copyfile(file2img,newdir);
           copyfile(file2hdr,newdir);

           % now prefix the copied file with the counter, and rename
           qq = counter(f);
           fnameIfin = strcat('file_',char(qq),'_',fnameI);
           fnameHfin = strcat('file_',char(qq),'_',fnameH);

           movefile(fnameI,fnameIfin);
           movefile(fnameH,fnameHfin);

        elseif c == 1 % some other file type

           % just copy the existing file
           copyfile(file,newdir);   

           % now prefix the copied file with the counter, and rename
           qq = counter(f);
           fnamefin = strcat('file_',char(qq),'_',fname);

           movefile(fname,fnamefin);

        end


    end

elseif typ == 2

    for f = 1:numf
        
        % get this file's target directory; it will already have filesep at
        % the end
        newdir = char(alltar(f,:));

        % get this file
        file = files(f,:);
        file = strcat(file); % to get rid of trailing spaces
        fend = numel(file);

        % get the location where th file name starts
        len = max(strfind(file,filesep));

        % do we have a ",1" tag?
           if strcmp(file(fend-1:fend),',1')
              file = file(1:fend-2);
              fend = fend - 2; % respecify this
           end

        % if dealing with .hdr/.img, prepare the .img
        
        suff = file(fend-3:fend); % the suffix of the file (assume .xxx)
        if strcmp(suff,'.img')
            %fname = file(len+1:fend-3)  % the actual file name after the last divider
            c = 0;
        else % any other format, including .nii, .xls, etc
            c = 1;    
        end   

        %fname = file(len+1:fend);  % the actual file name after the last divider

        if c == 0 % .hdr .img format

           file2 = file(1:fend-4); % get rid of ".img" 

           fileI = strcat(file2,'.img'); % full path of original img
           fileH = strcat(file2,'.hdr'); % full path of original header
           
           % now we either rename the file, or keep the existing name
           if changename == 1
               % rename it
               fileInew = strcat(newdir,newfilename,'.img');
               fileHnew = strcat(newdir,newfilename,'.hdr');
           elseif changename == 2
               % keep the original name
               fname = file(len+1:fend-4); % just the file name, without ".img"
               fileInew = strcat(newdir,fname,'.img');
               fileHnew = strcat(newdir,fname,'.hdr');
           end
           
           % and now copy to the new directory
           copyfile(fileI,fileInew);
           copyfile(fileH,fileHnew);

        elseif c == 1 % some other file type, including .nii

            % now we either rename the file, or keep the existing name
           if changename == 1
               % rename it
               filenew = strcat(newdir,newfilename,suff);
           elseif changename == 2
               % keep the original name
               fname = file(len+1:fend-4); % just the file name, without ".xxx"
               filenew = strcat(newdir,fname,suff); % add the suffix back on
           end
           
          % now copy the original file to its new name
           copyfile(file,filenew);  

        end


    end
    
    
end

fprintf(' Complete. \n\n');




