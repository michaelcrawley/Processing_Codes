% build the mpgread.dll MEX-file
%

mex -DWIN32 mpgread.c util.c video.c parseblk.c mvector.c decoders.c ...
				gdith.c jrevdct.c ordered2.c
            
            
