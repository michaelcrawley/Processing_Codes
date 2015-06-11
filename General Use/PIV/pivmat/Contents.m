% PIVMat Toolbox
% Version 3.02 25-June-2013
% F. Moisy
%
% Data import and display
%   loadvec           - Load vector/image fields
%   showf             - Display vector/scalar fields (still images or movies)
%   vec2scal          - Convert vector into scalar fields
%   imvectomovie      - Create a movie directly from files
%
% Field processing (both for vector and scalar fields)
%   filterf           - Spatial filter of a field 
%   interpf           - Interpolate missing data
%   bwfilterf         - Spatial Butterworth filter of a field 
%   medianf           - Median filter
%   averf             - Average of a series of fields 
%   spaverf           - Spatial average over X and/or Y of a field 
%   azaverf           - Azimuthally average field
%   phaseaverf        - Phase average fields
%   subaverf          - Substract the average of field 
%   smoothf           - Temporal smooth of fields 
%   resamplef         - Re-sample time series of fields
%   remapf            - Remap a field to a new grid
%   truncf            - Truncate a field 
%   extractf          - Extract a rectangle from a field 
%   resizef           - Resize a field
%   rotatef           - Rotate a field about the center 
%   shiftf            - Shift the axis of a field 
%   setoriginf        - Set the origin (0,0) of a field 
%   flipf             - Flip (mirror) a field
%   changefieldf      - Change the values of variables in a PIVMat structure
%   operf             - Other operations on fields 
%   addnoisef         - Add noise to a field
%   gradientf         - Gradient of a scalar field
%   nam               - Normalized Angular Momentum
%   matrixcoordf      - Matrix-coordinates of a point in a field
%   convert3dto2df    - Convert 3D vector fields to 2D vector fields
%
% Histograms, spectra and other statistics
%   statf             - Mean, standard deviation and rms of a field.
%   histf             - Histogram of a field.
%   histvec_disp      - Display histograms from vector fields.
%   histscal_disp     - Display histograms from scalar fields.
%   statvec_disp      - Display statistics from vector fields.
%   corrf             - Spatial correlation function of a scalar field.
%   tempcorrf         - Temporal correlation function of vector/scalar fields
%   vsf               - Structure functions of vector fields.
%   ssf               - Structure functions of scalar fields.
%   vsf_disp          - Display structure functions.
%   specf             - 1D power spectrum of vector/scalar fields.
%   spec2f            - 2D power spectrum of vector/scalar fields.
%   tempspecf         - Temporal power spectrum of vector/scalar fields.
%   jpdfscal          - Joint probability density function
%   jpdfscal_disp     - Display joint PDF
%   stresstensor      - Reynolds Stress tensor
%
% Advanced file operations
%   loadpivtxt        - Load a DaVis vector field exported in TXT mode
%   loadarrayvec      - Load a 2D array of vector fields
%   batchf            - Execute some functions for a series of files
%   vec2mat           - Convert DaVis files into MAT-Files
%   zerotonanfield    - Convert missing data representation to NaNs
%   nantozerofield    - Convert missing data representation to zeros
% 
% VEC/VC7 and SET Attribute handling
%   getattribute      - Get attributes from an IMX or VEC/VC7 file
%   readsetfile       - Read attributes from a .SET or .EXP file
%   getpivtime        - Acquisition time of a IMX/IM7 or VEC/VC7 files
%   getframedt        - Time intervals between the frames of an IMX/IM7 file
%
% Surface-Gradient Background-oriented Schlieren (FS-SS)
%   makebospattern    - Make a random dot pattern for SS or BOS
%   surfheight        - Surface height reconstruction.
%
% Miscellaneous files handling (from Fileseries v1.40)
%   cdw               - Change current working directory with wildcards
%   lsw               - List directory with wildcards
%   rdir              - Recursive list directory.
%   rdelete           - Delete files recursively.
%   rrmdir            - Delete directories recursively.
%   renamefile        - Rename a series of files.
%   renumberfile      - Re-number the indices of a series of files
%   getfilenum        - Get the index of a series of files.
%
% Miscellaneous
%   getsetname        - Get the name of the current set
%   getvar            - Get the value of the parameters in a string 'p1=v1_p2=v2_...'
%   randvec           - Random vector field.
%   vortex            - Vector field with a centered vortex
%   multivortex       - Vector fields of randomly distributed Burgers vortices

