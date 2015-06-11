function [ev rms fr]=averf(v,opt)
%AVERF  Average of vector/scalar fields
%   AF = AVERF(F) returns the average of the vector/scalar fields F.
%   AF is a field of the same type as F, whose elements are the average
%   of the elements of the fields F.
%
%   Depending on the nature of the fields F, the field AVERF(F) is usually
%   called "Ensemble average", "Phase average",  or "Time average".
%
%   By default, AVERF considers that the zero elements of F are erroneous,
%   and does not include them in the computations. If however you want to
%   force the zero elements to be included in the computations, specify
%   AVERF(F,'0').
%
%   [AF RMS FR] = AVERF(F,..) also returns the RMS field and the
%   fluctuation rate field of the vector/scalar fields F.
%      - RMS is a scalar field, whose elements are the rms of the elements
%        of the fields F.
%      - FR is a scalar field, whose elements are the local RMS normalized
%        by the local norm of AF.
%
%   If no output argument is specified, the average field AF is displayed.
%
%   Examples:
%      v = loadvec('*.vc7');
%      showf(averf(v));
%
%      [af rms fr]=averf(vec2scal(v,'uy'));
%      showf(fr,'clim',[0 1]);
%
%      averf *.vc7      % shows the average of the fields matching *.vc7
%
%   See also SMOOTHF, SPAVERF, SUBAVERF, AZAVERF, PHASEAVERF


%   F. Moisy, moisy_at_fast.u-psud.fr
%   Revision: 1.30,  Date: 2013/05/09
%   This function is part of the PIVMat Toolbox


% History:
% 2005/02/06: v1.00, first version.
% 2005/10/11: v1.10, averf (from ensavervec v1.00, 2005/10/06). also
%                    computes the rms field.
% 2006/03/10: v1.11, new history
% 2006/04/25: v1.12, bug fixed for scalar fields
% 2013/02/22: v1.20, works with 3D fields
% 2013/05/08: v1.30, bug fixed for computation of rms with 0 components
%                    (thanks Maciej Bujalski)

error(nargchk(1,2,nargin));
if (ischar(v) || iscellstr(v) || isnumeric(v))
    v=loadvec(v);
end

if nargin<=1, opt=''; end

ev=v(1);

vecmode=isfield(v(1),'vx');   % 1 for vector field, 0 for scalar field

nx=length(v(1).x);
ny=length(v(1).y);

nf = length(v); % number of fields

if vecmode
    ev.vx=zeros(nx,ny);
    ev.vy=zeros(nx,ny);
    nzx=zeros(nx,ny);
    nzy=zeros(nx,ny);
    if isfield(v(1),'vz')
        ev.vz=zeros(nx,ny);
        nzz=zeros(nx,ny);
    end
    for n=1:nf,
        ev.vx = ev.vx+v(n).vx;
        ev.vy = ev.vy+v(n).vy;
        nzx = nzx+logical(v(n).vx~=0); % counts the nonzero elements
        nzy = nzy+logical(v(n).vy~=0); % counts the nonzero elements
        if isfield(v(n),'vz')
            ev.vz = ev.vz+v(n).vz;
            nzz = nzz+logical(v(n).vz~=0);
        end
    end
    
    if strfind(opt,'0'),
        ev.vx=ev.vx/nf; % normalise with all the elements (including 0)
        ev.vy=ev.vy/nf;
        if isfield(v(1),'vz')
            ev.vz=ev.vz/nf;
        end
    else
        ev.vx = ev.vx./nzx;
        ev.vx(isnan(ev.vx)) = 0;
        ev.vy = ev.vy./nzy;
        ev.vy(isnan(ev.vy)) = 0;
        if isfield(v(1),'vz')
            ev.vz = ev.vz./nzz;
            ev.vz(isnan(ev.vz)) = 0;
        end
    end
    
else
    
    ev.w=zeros(nx,ny);
    nzw=zeros(nx,ny);
    for n=1:nf,
        nzw = nzw+logical(v(n).w~=0);
        ev.w = ev.w+double(v(n).w); % changed v1.12
    end
    if strfind(opt,'0')  % bug fixed here, v1.20
        ev.w = ev.w/nf;
    else
        ev.w = ev.w./nzw;
        ev.w(isnan(ev.w)) = 0;
    end
end

ev.name=['<' v(1).name '...' v(end).name '>'];

ev.history = {{v.history}' ['averf(ans, ''' opt ''')']}';

if nargout==0
    showf(ev);
    clear ev;
end

if nargout>1
    rms=v(1);
    if vecmode
        rms=rmfield(rms,{'vx','vy','unitvx','unitvy','namevx','namevy'});
        if isfield(v(1),'vz')
            rms=rmfield(rms,{'vz','unitvz','namevz'});
        end
        rms.unitw = v(1).unitvx;
        rms.namew = 'rms(u)';
    else
        rms.namew = ['rms(' v(1).namew ')'];
    end
    rms.w = zeros(nx,ny);
    nzw=zeros(nx,ny);
    for n=1:nf
        if vecmode
            if isfield(v(1),'vz')
                rms.w = rms.w + ((v(n).vx - ev.vx).*logical(v(n).vx)).^2 + ((v(n).vy - ev.vy).*logical(v(n).vy)).^2 + ((v(n).vz - ev.vz).*logical(v(n).vz)).^2  ; %FIXED
                nzw = nzw+logical((v(n).vx~=0)&(v(n).vy~=0)&(v(n).vz~=0));
            else
                rms.w = rms.w + ((v(n).vx - ev.vx).*logical(v(n).vx)).^2 + ((v(n).vy - ev.vy).*logical(v(n).vy)).^2 ; % FIXED
                nzw = nzw+logical((v(n).vx~=0)&(v(n).vy~=0));
            end
        else
            rms.w = rms.w + ((v(n).w - ev.w).*logical(v(n).w)).^2; % FIXED
            nzw=nzw+logical(v(n).w~=0);
        end
    end
    if strfind(opt,'0')
        rms.w = rms.w/nf;
    else
        rms.w = rms.w./nzw;
    end
    rms.w = rms.w.^0.5;
    rms.history = {{v.history}' 'rms(ans)'}';
end

if nargout>2
    fr=rms;
    if vecmode
        fr.namew = 'fluct-rate(u)';
    else
        fr.namew = ['fluct-rate(' v(1).namew ')'];
    end
    fr.unitw='';
    for i=1:nx
        for j=1:ny
            if vecmode
                if isfield(v(1),'vz')
                    if ((ev.vx(i,j)~=0) && (ev.vy(i,j)~=0)) && (ev.vz(i,j)~=0)
                        fr.w(i,j)=rms.w(i,j)/((ev.vx(i,j)^2+ev.vy(i,j)^2+ev.vz(i,j)^2)^(0.5));
                    end
                else
                    if ((ev.vx(i,j)~=0) && (ev.vy(i,j)~=0)),
                        fr.w(i,j)=rms.w(i,j)/((ev.vx(i,j)^2+ev.vy(i,j)^2)^(0.5));
                    end
                end
            else
                if ev.w(i,j)~=0,
                    fr.w(i,j)=rms.w(i,j)/abs(ev.w(i,j));
                end
            end
        end
    end
    fr.history = {{v.history}' 'fluct-rate(ans)'}';
end
