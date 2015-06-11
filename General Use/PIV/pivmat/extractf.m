function ef=extractf(f,inpcoord,opt)
%EXTRACTF  Extract a rectangular area from a vector/scalar field.
%   EF = EXTRACTF(F,[X1 Y1 X2 Y2]) extracts a rectangular area of
%   coordinates [X1 Y1 X2 Y2] from the original vector/scalar field(s) F.
%   By default, the coordinates are given in mesh units. Specify
%   EXTRACTF(F,[X1 Y1 X2 Y2],'phys') to give them in physical units (e.g.
%   in mm).
%
%   If no output argument, the result is displayed by SHOWF.
%
%   See also TRUNCF, ROTATEF, SHIFTF, FLIPF, REMAPF.


%   F. Moisy, moisy_at_fast.u-psud.fr
%   Revision: 1.40,  Date: 2013/02/22
%   This function is part of the PIVMat Toolbox


% History:
% 2005/02/01: v1.00, first version.
% 2005/02/23: v1.01, cosmetics.
% 2005/07/26: v1.02, added history.
% 2005/09/06: v1.03, check arg added.
% 2005/09/15: v1.04, no recursive call.
% 2005/10/11: v1.10, now called extractf, from extractvec 1.04 (2005/09/15)
% 2006/07/21: v1.20, now works with physical units or mesh units.
% 2010/05/04: v1.30, accepts coordinates outside the physical bounds
% 2013/02/22: v1.40, works with 3D fields

error(nargchk(2,3,nargin));

if (ischar(f) || iscellstr(f) || isnumeric(f))
    f=loadvec(f);
end 

if nargin<3,
    opt='mesh';  % center given in physical units by default
end

if strncmpi(opt,'phys',1),    
    % converts physical units into mesh units
    coord(1) = 1 + floor((inpcoord(1) - f(1).x(1))/abs(f(1).x(2)-f(1).x(1)));
    coord(3) = 1 + ceil((inpcoord(3) - f(1).x(1))/abs(f(1).x(2)-f(1).x(1)));
    
    coord(2) = 1 + floor((inpcoord(2) - f(1).y(1))/abs(f(1).y(2)-f(1).y(1)));
    coord(4) = 1 + ceil((inpcoord(4) - f(1).y(1))/abs(f(1).y(2)-f(1).y(1)));
else
    coord = inpcoord;
end

nx=length(f(1).x);
ny=length(f(1).y);

coord(1) = max([coord(1) 1]);
coord(3) = min([coord(3) nx]);
coord(2) = max([coord(2) 1]);
coord(4) = min([coord(4) ny]);
    
vecmode=isfield(f(1),'vx');   % 1 for vector field, 0 for scalar field

ef=f;

for i=1:length(f)
    if vecmode
        ef(i).vx = f(i).vx(coord(1):coord(3),coord(2):coord(4));
        ef(i).vy = f(i).vy(coord(1):coord(3),coord(2):coord(4));
        if isfield(f(i),'vz')
            ef(i).vz = f(i).vz(coord(1):coord(3),coord(2):coord(4));
        end
    else
        ef(i).w = f(i).w(coord(1):coord(3),coord(2):coord(4));
    end
    ef(i).x = f(i).x(coord(1):coord(3));
    ef(i).y = f(i).y(coord(2):coord(4));
    ef(i).history = {f(i).history{:} ['extractf(ans, [' num2str(inpcoord) '], ''' opt ''')']}';
end

if nargout==0
    showf(ef);
    clear ef
end
