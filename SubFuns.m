function SubFuns(mFile)
%SUBFUNS lists all function declaration lines in specified MFILE.
%   SUBFUNS(MFILE) displays list to Command Window.
%
%Example:
%   SubFuns SubFuns
%
%Note:
%   This utility uses both subfunctions and a nested function, both of
%   which are supported.
%
%See also FUNCTION.
%
% Copyright 2004-2010 The MathWorks, Inc.
%
% http://www.mathworks.com/matlabcentral/fileexchange/5563-subfuns#

%{
Introducing SubFuns, a command-line utility to list all function declarations within a specified m-file.

Contents

MATLAB functions
Different function declarations
Problem: How many functions inside an m-file?
Solution: custom command-line utility
Example
MATLAB functions

Functions are an integral part of MATLAB. They provide higher levels of functionality and abstraction. For example, the standard deviation of a matrix is the square root of the variance, where the variance is the sum of the squares of the elements divided by the total number of elements. So, in MATLAB the std function calls the var function, which in turn calls the sum function.

This is just one of many examples where MATLAB makes extensive use of functions. Another is toolboxes available from The MathWorks and third parties, which are largely just collections of trusted, proven, standard functions.

In addition, many MATLAB users write their own custom functions and toolboxes.

Different function declarations

MATLAB 7 allows three types of function declarations within an m-file. To learn more, see the associated online documentation.

Primary M-File Function http://www.mathworks.com/access/helpdesk/help/techdoc/matlab_prog/ch_func8.shtml
Subfunctions http://www.mathworks.com/access/helpdesk/help/techdoc/matlab_prog/ch_fun14.shtml
Nested Functions http://www.mathworks.com/access/helpdesk/help/techdoc/matlab_prog/ch_func9.shtml
Problem: How many functions inside an m-file?

With all this programming flexibility, it's not uncommon for an m-file to contain numerous function declarations. The MATLAB Editor provides a pull-down tool to Show Functions within an open m-file.

Suppose you want to see function decalarations inside an m-file without editing. Perhaps you need to debug GUI callbacks, or understand someone else's function. You may not remember everything in a long list of functions after the Show Functions pull-down disappears. In these situations it would be helpful if the function list would hang around for reference.

Solution: custom command-line utility

To solve this problem, I created a command-line utility function, called SubFuns. Like all good MATLAB programs, SubFuns includes Help text.

help SubFuns
 SUBFUNS lists all function declaration lines in specified MFILE.
    SUBFUNS(MFILE) displays list to Command Window.
 
 Example:
    SubFuns SubFuns
 
 Note:
    This utility uses both subfunctions and a nested function, both of
    which are supported.



Example

SubFuns is a MATLAB function that contains all three kinds of function declarations. To see these functions we can run SubFuns on itself.

SubFuns SubFuns
   Line:   Function:
   1       SubFuns
   77      FunctionName
   94      NoLeadBlanks
   152     OpenLink


Published with MATLAB® 

%}



error(nargchk(1,1,nargin))
error(nargoutchk(0,0,nargout))

%make sure m-file exists (on path)
if ~exist(mFile)
  error(['Unable to find m-file ' mFile])
end

%open m-file (read access)
fid = fopen(which(mFile));

%find all lines starting with 'function' signature
subFunList = {};
lineNum = 0;
while ~feof(fid) %check one line at a time
  lineNum = lineNum+1; %advance line pointer
  try %read line
    txt = fgetl(fid);
  catch %file read error
    disp(['Tried reading file ' mFile])
    error(sprintf('Unable to read line %d',lineNum))
  end
  name = FunctionName(txt);
  if ~isempty(name) %found one, add to list
    subFunList(end+1,1:2) = {lineNum,name};
  end
end

%close file when done
fclose(fid);

%display results
if ~isempty(subFunList)
  %honor loose preference
  if strcmp('loose',get(0,'FormatSpacing'))
    fprintf('\n');
  end

  %table header
  fprintf('   Line:   Function:\n');

  %row/col entries
  for i=1:size(subFunList,1)
    lineNum = subFunList{i,1};
    %fprintf('   %d',lineNum);
    html = OpenLink(mFile,lineNum);
    fprintf('   %s',html);
    fprintf(blanks(7-floor(log10(lineNum))));
    name = subFunList{i,2};
    fprintf('%s\n',name);
  end

  %extra blank line after list before prompt
  fprintf('\n');
end %if (list empty or not)

end %SubFuns (main)


%sub-------------------------------------------------------------------
function name = FunctionName(txt)
%FUNCTIONNAME returns name of subfunction.
%   Returns empty if line does not contain valid function name.
%
%Usage:
%   name = FunctionName(txt)

%TODO List:
%   1. double check that algorithm worked correctly

%Note: MATLAB 7 includes a new TODO/FIXME tool to help maintain your code.
%   It's one of several new Directory Reports available from the Current
%   Directory Browser on the MATLAB Desktop.  You can also access the
%   TODO/FIXME Report tool from the menu bar when the Current Directory is
%   selected using View, Directory Reports.

error(nargchk(1,1,nargin))
error(nargoutchk(1,1,nargout))

%remove leading blanks
txt = NoLeadBlanks(txt);

  %nested---------------------------------------%
  function txt = NoLeadBlanks(txt)              %
  %NOLEADBLANK removes leading spaces.          %
    done = false;                               %
    while ~done                                 %
      if isempty(txt)                           %
        return                                  %
      elseif txt(1)==' ' %leading blank         %
        txt(1) = []; %remove                    %
      else %done (nonblank first character)     %
        done = true;                            %
      end                                       %
    end %while (done or not)                    %
  end %NoLeadBlanks (nested function)-----------%

  %{
    Notes:
     (1) This use of a nested function was contrived. There is no advantage
         over a subfunction in terms of variable scope/sharing. Its value 
         here was just to provide a self referential example to show a 
         nested fuction.
     (2) This is also an example of a block comment (new for MATLAB 7).
  %}
  
if length(txt)<8 %line too short
  name = []; %return empty (no valid function here)
  return;
end

if strcmp('function',txt(1:8)) %function line
  %strip off function keyword at beginning
  txt(1:8) = [];

  %ignore outputs
  if strfind(txt,'=') %outputs defined
    loc = strfind(txt,'=');     %LHS before =
    txt(1:loc) = [];            %ignore outputs
  end

  %remove leading blanks before name
  txt = NoLeadBlanks(txt);

  %ignore EOL comments
  if strfind(txt,'%')
    loc = strfind(txt,'%');     %comment from % on
    txt(loc:end) = [];          %ignore inputs
  end
  
  %ignore inputs
  if strfind(txt,'(')
    loc = strfind(txt,'(');     %inputs after (
    txt(loc:end) = [];          %ignore inputs
  end

  %remove trailing blanks after function name
  while txt(end)==' '
    txt(end) = [];
  end

  %for legal function names, only thing left must be function name
  name  = txt;
else %not function line
  name = []; %return empty
end

end %FunctionName (subfunction)


%sub-------------------------------------------------------------------
function html = OpenLink(mFile,lineNum)
%OPENLINK generates hyperlink to open MFILE to LINENUM.
%   HTML = OPENLINK(MFILE,LINENUM) returns HTML code to generate hyperlink.
%   Displays line number. Link invokes command to open specified MFILE to
%   LINENUM.
%
%Example:
%   html = OpenLink(mFile,lineNum);

error(nargchk(2,2,nargin))
error(nargoutchk(1,1,nargout))

%make m-file location explicit
mFile = which(mFile);

%Java function
javaFcn = 'com.mathworks.mlservices.MLEditorServices.openDocumentToLine';

%calling syntax
cmd = sprintf('matlab: %s(''%s'',%d,1,1)',javaFcn,mFile,lineNum);

%hyperlink HTML script
html = sprintf('<a href="%s">%d</a>',cmd,lineNum);

end %OpenLink (subfunction)
