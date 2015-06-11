function res = checkreadimxversion(a,ver)
%CHECKIMXVERSION  Check the version of the ReadIMX DLL
%   CHECKIMXVERSION(A,VER), where A is a structure obtained from READIMX,
%   returns TRUE if the DLL Version is larger or equal to VER.
%
%   F. Moisy
%   Revision: 1.00,  Date: 2008/04/04


res = false;
if isfield(a,'DLLVersion')
    str = a.DLLVersion;
    pos1 = findstr(str,'Revision: ');
    pos2 = findstr(str,'$'); pos2=pos2(end);
    curver = str2double(str((pos1+10):(pos2-1)));
    if ~isempty(curver)
        res = (curver>=ver);
    end
end

    