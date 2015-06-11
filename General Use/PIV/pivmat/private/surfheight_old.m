function hresult = surfheight(v,h0,n,mode,varargin)
%SURFHEIGHT Surface height
%   H = SURFHEIGHT(R, H0, N) computes the surface height H, given the
%   displacement field R, the mean height H0 and the refraction index N
%   (N=1.33 by default). H0 is expressed in the same units as R (e.g., mm).
%   The displacement field is obtained from the correlation of a reference
%   image (flat surface) and a distorted image (wavy surface).
%
%   H = SURFHEIGHT(R, H0, N, MODE) specifies the integration mode:
%    1  : WS      weak slope and arbitrary amplitude (by default)
%    2  : WS-WA   weak slope and weak amplitude
%    3  : WA      arbitrary slope and weak amplitude
%
%   H = SURFHEIGHT(..., 'submean') first subtracts the mean displacement
%   along each direction. This is useful if the incidence angle of the
%   camera is not exactly 0.
%
%   By default, SURFHEIGHT produces a height field such that MEAN(H)=H0.
%   Specify H = SURFHEIGHT(..., 'nosetzero') if you do not want to apply
%   this constraint - in this case, the point (1,1) of the height field
%   is equal to H0, but MEAN(H) is arbitrary.
%
%   H = SURFHEIGHT(..., 'verbose') displays the computation progress.
%
%   SURFHEIGHT is based on INTGRAD2.M by J. D'Errico.
%
%   Example:
%
%     r = loadvec('*.vc7');
%     h = surfheight(r, 10);
%     showscal(h);
%
%   See also SHOWVEC, SHOWSCAL, VEC2SCAL, GRADIENTF.

%   F. Moisy, moisy_at_fast.u-psud.fr
%   Revision: 1.31,  Date: 2007/07/19
%   This function is part of the PIVMat Toolbox

% History:
% 2006/07/01: v1.00, first version, Joran Rolland
% 2007/04/10: v1.10, help, arg check. Inclusion in PIVMat
% 2007/04/15: v1.11, new option 'mode'
% 2007/05/21: v1.12, new mode 3: arbitrary slope and small amplitude
% 2007/06/15: v1.20, new options 'submean' and 'setzero'
% 2007/07/16: v1.30, new modes 4-5-6 (camera above)
% 2007/07/19: v1.31, new option 'verb'; new submean


error(nargchk(2,inf,nargin));

if (ischar(v) || iscellstr(v) || isnumeric(v))
    v=loadvec(v);
end

if nargin<3
    n = 1.33;   % water index
end

if nargin<4
    mode=1; % default mode = small slope and arbitrary amplitude
end

if mode<=3
    alpha = 1-1/n;
else
    alpha = n-1;
end

scale = abs(v(1).x(1) - v(1).x(2));
[nx,ny] = size(v(1).vx);

hresult = vec2scal(v,'norm'); % final result : h

for i=1:length(v)
   

    if any(strncmpi(varargin,'submean',4))   % new v1.20
        % subtract the displacement averaged line-by-line:
        %dx = v(i).vx - ones(nx,1)*mean(v(i).vx,1);
        %dy = v(i).vy - mean(v(i).vy,2)*ones(1,ny);
              
        % subtract the displacement averaged on the whole image:
        dx = v(i).vx - mean(v(i).vx(:));
        dy = v(i).vy - mean(v(i).vy(:));
        
        % subtract the displacement averaged on part of the image:
        %dx = v(i).vx - mean(mean(v(i).vx(1:round(end/2),:)));
        %dy = v(i).vy - mean(mean(v(i).vy(1:round(end/2),:)));
    else
        dx = v(i).vx;
        dy = v(i).vy;
    end

    switch mode
        case {1,4}  % small slope and arbitrary amplitude
            dh2dx = -(2/alpha)*dx;
            dh2dy = -(2/alpha)*dy;
            h = intgrad2(dh2dy,dh2dx,scale,scale,h0^2).^(1/2);
            
        case {2,5}  % small slope and small amplitude
            dhdx = -(1/(alpha*h0))*dx;
            dhdy = -(1/(alpha*h0))*dy;
            h = intgrad2(dhdy,dhdx,scale,scale,h0);
            
        case 3  % arbitrary slope and small amplitude
            dhdx = -invphi(dx/h0, n);
            dhdy = -invphi(dy/h0, n);
            h = intgrad2(dhdy,dhdx,scale,scale,h0);
            
        case 6  % arbitrary slope and small amplitude - from below
            %dhdx = -invphi(dx/h0, n);
            %dhdy = -invphi(dy/h0, n);
            %h = intgrad2(dhdy,dhdx,scale,scale,h0);
    end

    if ~any(strncmpi(varargin,'nosetzero',4))   % new v1.20
        % re-set the mean height to h0:
        h = h - mean(h(:)) + h0;
    end
    
    hresult(i).w = h;
    hresult(i).namew = 'h';
    hresult(i).unitw = v(i).unitx;
    hresult(i).history = {v(i).history{:} ['surfheight(ans, ' num2str(h0) ', ' num2str(n) ', ' num2str(mode) ', ' varargin{:} ')']}';

    if any(strncmpi(varargin,'verbose',4))
        disp([num2str(i) '/' num2str(length(v)) ' done']);
    end
end

if nargout==0
    showscal(hresult);
    clear hresult
end
