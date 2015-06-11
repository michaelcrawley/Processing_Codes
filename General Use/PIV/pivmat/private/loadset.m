function vv=loadset(setname,filename,opt)
%LOADSET  Loads set(s) of vector/image fields (obsolete)
%
%      WARNING:  LOADSET is present in the private folder of the PIVMat
%                toolbox for compatibility issues. Use LOADVEC instead.
%
%   F = LOADSET loads all the vector/image fields of the current directory
%   into the structure array F. This is essentially equivalent to
%   F = LOADVEC({'*.vec','*.vc7','*.imx','*.im7'}), except that LOADSET
%   makes use of a MAT-file called 'set.mat' to save time (see below).
%   See LOADVEC for the content of the structure array F.   
%
%   F = LOADSET(DIR) loads all the vector/image fields from the specified
%   directory DIR in the structure array F. If DIR is a cell array
%   (e.g. {'dir1','dir2'}), then loads all the files from each directory
%   and concatenates them in a single structure array F (see LOADARRAYVEC
%   to load the files in a 2D structure array). Wildcards (*) may be used:
%   F = LOADSET('set*') loads all the files in all the directories matching
%   'set*'. Brackets ([]) are also accepted (see EXPANDSTR for details).
%   For example, F = LOADSET('set[1:3,1]') loads all the vector/scalar
%   fields from the directories 'set1', 'set2', 'set3'. DIR can be a cell
%   array with a combination of wildcards and brackets.
%   
%   F = LOADSET(DIR, FILE) loads only the file(s) FILE in the directory(s)
%   DIR. This is useful when DIR and/or FILE are cell arrays and/or contain
%   wildcards (*) and/or brackets (see EXPANDSTR). This is equivalent to
%   F = LOADVEC([DIR '/' FILE]).
%
%   F = LOADSET(DIR, NUM) loads only the file number NUM in the directory(s)
%   DIR (works in alphanumeric order, only for VEC/VC7 and IMX/IM7 files).
%   NUM may be a simple number or any valid MATLAB vector (e.g., 1:10,
%   [1 10], 4:-1:1, etc.)
%
%   LOADSET(...) without output argument is a shortcut for
%   SHOWVEC(LOADSET(...)) or SHOWSCAL(LOADSET(...))
%
%   Examples:
%
%     V = LOADSET('myset') loads all the VEC/VC7 files in the directory
%     'myset' (loads the 'set.mat' file if it exists; see below).
%
%     SHOWVEC(LOADSET) displays all the files from the current directory.
%
%     V = LOADSET('set*','B00001.vec') loads the files 'B00001.vec'
%     contained in each directory 'set*'.
%
%     V = LOADSET('set*','*.vec') loads all the VEC files in each directory
%     'set*'.
%
%     V = LOADSET({'set1','set2'},{'file1.vec','B*.vec'}) loads 'file1.vec'
%     and all the VEC files 'B*' in the two directories 'set1' and 'set2'.
%
%     V = LOADSET('set*','B[5:10]*') loads 'B00005.vec' to 'B00010.vec' in
%     each directory 'set*' (see LOADVEC and EXPANDSTR for details).
%
%     V = LOADSET('set[200:5:400,2]*','B01.vec') loads the files 'B01.vec'
%     contained in each directory 'set200*','set205*','set210*' etc.
%
%     V = LOADSET('set*',1:10) loads the 10 first files from each directory
%     'set*'.
%
%   The first time LOADSET is used to load all the vector/scalar fields in a
%   directory (i.e. without the argument FILE), a MAT-file called 'set.mat'
%   containing the structure array F for all the files is created. The next
%   time LOADSET is used, this file 'set.mat' is loaded instead of each
%   original file, to save time (this may be 2 to 5 times faster). Once
%   this file 'set.mat' has been created, the original files are not
%   necessary any more.
%
%   F = LOADSET(...,'overwrite') loads the original vector/scalar fields
%   even if a 'set.mat' file already exists, and saves (overwrites) a new
%   'set.mat'.
%
%   F = LOADSET(...,'nosave') loads the original vector/scalar files even
%   if a 'set.mat' file already exists, and does not save the 'set.mat'
%   file.
%
%   See also LOADVEC, VEC2MAT, LOADARRAYVEC, BATCHF, READSETFILE, SHOWVEC, RDIR.


%   F. Moisy, moisy_at_fast.u-psud.fr
%   Revision: 1.50,  Date: 2010/05/03
%   This function is part of the PIVMat Toolbox


% History:
% 2004/04/21: v1.00, first version.
% 2005/02/07: v1.01, added the 'setname' field.
% 2005/02/16: v1.02, now works with arbitrary file numbers
% 2005/02/22: v1.03, cosmetics
% 2005/04/22: v1.10, allows to load multiple sets.
% 2005/04/28: v1.11, added the 'history' field.
% 2005/06/07: v1.20, Compatible VEC/VC7 (DaVis 6/7).
% 2005/09/01: v1.21, cosmetics.
% 2005/09/06: v1.22, use cell arrays instead of string arrays for multiple sets.
% 2005/09/09: v1.30, wildcards allowed, and option 'overwrite' added.
% 2005/09/28: v1.40, allows to load several files in each set.
%                    option 'nosave' added. loadset without output arg.
%                    added
% 2005/09/29: v1.41, filenames are automatically expanded if needed.
% 2005/09/30: v1.42, setnames are automatically expanded too.
% 2005/10/06: v1.43, cosmetics.
% 2005/10/19: v1.44, help text changed.
% 2005/11/18: v1.45, also accept file numbers (only help txt changed)
% 2005/12/20: v1.46, help text changed.
% 2006/03/03: v1.47, history is now given by loadvec's history
% 2006/03/10: v1.48, imx/im7 files accepted
% 2006/04/26: v1.49, bug fixed in the input argument processing
% 2010/05/03: v1.50, new warning message



warning('PIVMat:loadset:obsolete',...
    'LOADSET will be obsolete in future versions of PIVMAT. Use LOADVEC instead.');

if ~exist('setname','var'), setname=''; end
if ~exist('filename','var'), filename=''; end
if ~exist('opt','var'), opt=''; end

% reorder the arguments if setname and/or filename is absent:
if (strcmpi(setname,'overwrite')) || (strcmpi(setname,'nosave'))
    opt=setname; setname=''; filename=''; % bug fixed v1.49
end

if (strcmpi(filename,'overwrite')) || (strcmpi(filename,'nosave'))
    opt=filename; filename='';
end

%expand the setname if necessary (v1.42):
if ~iscell(setname)
    setname=expandstr(setname);
end

if isempty(setname)
    %if no setname specified, loads the current set:
    v=loadcurset(opt);  %(this internal function is defined below)
else
    if (~iscell(setname)) && any(strfind(setname,'*'))
        % if setname is not a cell array, check for wildcards:
        listsetname=dir(setname);
        % setname is now a cell array containing all the
        % sets (directories) matching the wildcards:
        setname={listsetname.name};
    end
    if ~isempty(setname)
        if iscell(setname),
            %if a cell array is specified, calls loadset
            %for each element of the cell and concatenates the sets:
            v=[];
            for i=1:length(setname);
                v=[v loadset(setname{i},filename,opt)];
            end
        else
            % for a single setname, loads the set
            % first removes the '.set' suffix if present
            if length(setname)>=4,
                if isequal(setname((length(setname)-3):length(setname)),'.set'),
                    setname=setname(1:(length(setname)-4));
                end
            end
            if exist(setname,'dir')
                cd(setname);
                if ~isempty(filename), % new v1.40
                    % loads only the file(s) filename in this directory
                    v=loadvec(filename);
                else       
                    % loads all the file(s) in this directory
                    v=loadcurset(opt);   %(defined below)
                end
                cd ..
            else
                error(['Directory ''' setname ''' not found.']);
            end
        end
    end
end


if nargout==0
    showf(v);
else
    vv=v;
end

% -----------------------------------------------------------------------

function v=loadcurset(opt)
%  This internal subfunction simply loads the set of the current directory,
%  and creates the file 'set.mat' if it does not already exist.
if nargin==0
    opt='';
end
if (exist('set.mat','file')) && (~strcmpi(opt,'overwrite'))
    load set.mat;  % the set is stored in the variable 'v' of this file.
else
    vecfiles=[dir('*.vec') dir('*.vc7') dir('*.imx') dir('*.im7')];
    for i=1:length(vecfiles)
        v(i)=loadvec(vecfiles(i).name); % changes Feb 16, 2005
        v(i).setname=getsetname; % Feb 07, 2005 FM ; bug fixed June 7, 2005.
    end
    if ~strcmpi(opt,'nosave') % Sept 28, 2005, v1.31
        save('set.mat','v');
    end
end
