function mov = imvectomovie(filename, varargin)
%IMVECTOMOVIE  Convert a series of IM7 or VC7 files into a movie
%   MOV = IM7TOMOVIE(FILENAME, 'options', ...) converts the series of DaVis
%   files matching FILENAME into a movie MOV, that can be saved as an AVI
%   file using the command MOVIE2AVI.  All data formats supported by
%   LOADVEC are allowed (e.g., VC7, IM7 etc.) FILENAME is a string which
%   can contain wildcards (*) and brackets for file number enumeration
%   (see LOADVEC for all available syntax).
%
%   The Image Processing toolbox is required for this function.
%
%   This command is essentially similar to MOV = SHOWF(FILENAME,...).
%   The essential difference is that here the whole files are NOT stored
%   into a PIVMat structure, but are loaded one by one. This allows the
%   user to generate very large movies without 'out of memory' problems.
%
%   All the options of SHOWF (e.g., 'clim', 'cmap', 'surf', 'title' etc.)
%   are available with IMVECTOMOVIE.
%
%   Example
%      mov = imvectomovie('B*.im7', 'clim', 2000, 'verbose');
%      movie2avi(mov, 'mymovie.avi', 'fps', 10);
%
%   See also MOVIE2AVI, LOADVEC, SHOWF

%   F. Moisy, moisy_at_fast.u-psud.fr
%   Revision: 1.00,  Date: 2011/05/01
%   This function is part of the PIVMat Toolbox

% History:
% 2011/05/01: v1.00, first version.

file = rdir(filename,'fileonly'); % resolve wildcards
n = length(file);
for i=1:n
    if any(strncmpi(varargin,'verbose',4))
        disp(['Loading file #' num2str(i) '/' num2str(n) ': ' file{i}]);
    end    
    mov(i) = showf(loadvec(file{i}), varargin{:});
end
