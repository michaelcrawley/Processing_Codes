function [status,output] = Archiver(mode,archive,varargin)

    delchk = 0; %initialize value
    switch lower(mode)
        case {'a','archive'}
            cmd = ' a '; %specify command to 7za.exe
            
            %check whether user wants files deleted after archiving
            tmp = strcmpi(varargin,'-deletefiles');
            delchk = any(tmp);
            if delchk %remove this switch as it is a matlab switch, not dos switch
               varargin = varargin(~tmp);
            end
            
            %specify file list for compression
            flist = varargin{1}; 
            if ~iscell(flist), flist = {flist}; end
            arg = sprintf(' "%s" ',flist{:});
                        
            if length(varargin) > 1 %check options
                options = varargin(2:end); 
                if ~iscell(options), options = {options}; end
                
                %specify compression level
                cmpstr = {'x0','x1','x3','x5','x7','x9'};
                levelchk = cellfun(@(str) any(strcmpi(str,cmpstr)),options);
                if any(levelchk)
                    level = [' -m' options{levelchk},' '];
                    options = options(~levelchk);
                else
                    level = ' -mx9 ';
                end
                
                %specify compression algorithm - 7za appears to be non-case
                %dependent
                cmpstr = {'Copy','Deflate','Deflate64','BZip2','LZMA','LZMA2','PPMd','Delta','BCJ','BCJ2'};
                algchk = cellfun(@(str) any(strcmpi(str,cmpstr)),options);
                if any(algchk)
                    alg = [' -m0=' options{algchk},' '];
                    options = options(~algchk);
                end  
            else
                level = ' -mx9 ';
                options = {};
            end
            
            if ~exist('alg','var') %set default compression if not specified by user
                [~,~,ext] = fileparts(archive);
                switch ext
                    case '.7z'
                        alg = ' -m0=LZMA2 ';
                    case '.zip'
                        alg = ' -mm=LZMA ';
                    otherwise
                        alg = '';                        
                end                        
            end
            
            stdswt = ' -ssw -mtc=on'; %enables compression of open files, stores file creation timestamps
            
            %specify additional user defined switches
            if ~isempty(options)
               misc = sprintf(' %s ',options{:});
            else
                misc = '';
            end
            
            swt = [stdswt,level,alg,misc];
        case {'l','list'}
            cmd = ' l ';
            arg = '';
            
            %specify additional user defined switches
            if ~isempty(varargin)
               swt = sprintf(' "%s" ',varargin{:});
            else
                swt = '';
            end           
        case {'x','extract'}
            cmd = ' x ';
            arg = '';
            
            if ~isempty(varargin) %sets output directory
                dir_out = [' -o"',varargin{1},'" '];
            else
               dir_out = [' -o"',pwd,'" ']; 
            end
            
            stdswt = ' -aou -y '; %sets default to auto rename extracted files, assumes yes to all queries
            swt = [stdswt,dir_out];
    end
    
    [status,output] = dos(['7za',cmd,swt,archive,arg]);
     
    if delchk && status==0 %delete files if user desires and if they were properly archived
        type = cellfun(@(var) exist(var),flist);
        files = type==2;
        folders = type==7;
        cellfun(@(file) delete(file),flist(files));
        cellfun(@(folder) rmdir(folder,'s'),flist(folders));
    end
end

