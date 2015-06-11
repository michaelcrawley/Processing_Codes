function fv=bwfilterf(v,filtsize,order,varargin)
%BWFILTERF  Butterworth filter for a vector/scalar field
%   FF = BWFILTERF(F,FSIZE,ORDER) applies a lowpass Butterworth filter
%   to the vector/scalar field F  with a cutoff size FSIZE (in grid units)
%   and order ORDER. Dimensions of FSIZE and ORDER must match the dimension
%   of the field F. If FSIZE and/or ORDER are scalar, the same cutoff and
%   order is applied to each field of F. BWFILTERF first Fourier-transforms
%   the field(s), applies a low/high-pass transfer function and inverse
%   Fourier transforms. Typical values for FSIZE are around 1-10, and typical
%   values for ORDER are in the range 2-10 (positive for a low-pass filter,
%   negative for a high-pass filter). Sharp filters (ie large values of
%   FSIZE) may produce oscillations.
%
%   FF = BWFILTERF(F,FSIZE,ORDER,OPT), where OPT is a string that may
%   contain one or several substrings:
%     'low', 'high':   specifies a lowpass (by default) or highpass filter
%     'trunc':         truncates the borders of width FSIZE, which are
%                      affected by the filtering.
%
%   The X and Y dimensions of the fields must be even. If one of the
%   dimension is odd, the last column/row is discarded.
%
%   If no output argument, the result is displayed by SHOWF.
%
%   Note 1: A highpass filter of order ORDER is equivalent to a lowpass
%   filter of order -ORDER.
%
%   Note 2: If there are missing data in the field, it is better to first
%   interpolate the data. See INTERPF.
%
%   Example:
%      v = loadvec('*.vc7');
%      showf(bwfilterf(v,3,8));
%
%   See also FILTERF, INTERPF, ADDNOISEF, TRUNCF, EXTRACTF.


%   F. Moisy, moisy_at_fast.u-psud.fr
%   Revision: 1.40,  Date: 2013/02/22
%   This function is part of the PIVMat Toolbox


% History:
% 2005/02/23: v1.00, first version.
% 2005/05/25: v1.01, bug fixed.
% 2005/09/03: v1.02, arg check added.
% 2005/10/06: v1.03, highpass filter added.
% 2005/10/09: v1.10, now called bwfilterf (from bwfiltervec 1.03, 2005/10/06)
% 2006/03/10: v1.11, new history
% 2006/04/25: v1.12, bug fixed with trunc/truncf
% 2008/10/08: v1.20, now accepts non-square fields. Vectorized version.
% 2009/03/17: v1.30, FSIZE and ORDER can be arrays.
% 2013/02/22: v1.40, works with 3D fields


error(nargchk(3,4,nargin));
if (ischar(v) || iscellstr(v)), v=loadvec(v); end

if any(strncmpi(varargin,'high',1))
    order = -order;   % for a high-pass filter
end

if numel(filtsize)>1
    if numel(filtsize)~=numel(v)
        error('Dimensions of vector/scalar fields F and filter sizes must be equal.');
    end
else
    filtsize = filtsize * ones(1,numel(v));
end

if numel(order)>1
    if numel(order)~=numel(v)
        error('Dimensions of vector/scalar fields F and order must be equal.');
    end
else
    order = order * ones(1,numel(v));
end

nx=length(v(1).x);
ny=length(v(1).y);

% if size is even, truncate
if mod(nx,2)==1
    nx = nx-1;
    v = extractf(v, [1 1 nx ny]);
end

if mod(ny,2)==1
    ny = ny-1;
    v = extractf(v, [1 1 nx ny]);
end

n=min(nx,ny);  % note sure of this for rectangular fields...
kx0=nx/2+1; % location of the zero mode.
ky0=ny/2+1;

% computes the transfer function (spectral filter):
[ky,kx] = meshgrid(1:ny,1:nx);
k = sqrt((kx-kx0).^2+(ky-ky0).^2);

fv=v;

vecmode=isfield(v(1),'vx');   % 1 for vector field, 0 for scalar field

for i=1:numel(v)  % do NOT use parfor here!
    if filtsize(i)~=0
        T = 1./(1+(k/(n/filtsize(i))).^(order(i)/2));
        if vecmode,
            spx=fftshift(fft2(v(i).vx));
            spy=fftshift(fft2(v(i).vy));
            fv(i).vx = real(ifft2(ifftshift(spx.*T)));
            fv(i).vy = real(ifft2(ifftshift(spy.*T)));
            if isfield(v(i),'vz')
                spz=fftshift(fft2(v(i).vz));
                fv(i).vz = real(ifft2(ifftshift(spz.*T)));
            end
        else
            sp=fftshift(fft2(v(i).w));
            fv(i).w = real(ifft2(ifftshift(sp.*T)));
        end
        if any(strncmpi(varargin,'trunc',1))
            fv(i)=truncf(fv(i),floor(filtsize(i)));
        end
        fv(i).history = {v(i).history{:} ['bwfilterf(ans, ' num2str(filtsize(i)) ', ' num2str(order(i)) ', ''' varargin{:} ''')']}';
    end
end

if nargout==0
    showf(fv);
    clear fv
end
