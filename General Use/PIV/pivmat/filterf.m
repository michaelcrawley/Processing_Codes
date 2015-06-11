function fv=filterf(v,filtsize,method,varargin)
%FILTERF  Apply a lowpass filter to a vector/scalar field.
%   FF = FILTERF(F,FSIZE) applies a Gaussian filter of size FSIZE to the
%   vector/scalar field(s) F. If FSIZE is an array, its dimension must
%   match the dimension of F. If FSIZE is a scalar, the same filter size is
%   applied to each field of F. FSIZE=1 is taken by default if not specified.
%
%   FF = FILTERF(F,FSIZE,METHOD) specifies the filter:
%         'flat':    flat (or top-hat) matrix, ones(FSIZE,FSIZE)
%                    (FSIZE must be an even integer)
%         'gauss':   gaussian (by default)
%         'igauss':  derivative of the integrated gaussian (minimized
%                    discretisation effects for small FSIZE).
%
%   Note : If there are missing data in the field, it is better to first
%   interpolate the data. See INTERPF.
%
%   The size of the filtered field is smaller than the original field, to
%   avoid boundary effects (the convolution is done by CONV2 with the
%   option 'valid'). If you prefer to keep the whole field, use
%   FILTERF(F,FSIZE,METHOD,'same')
%
%   If no output argument, the result is displayed by SHOWF.
%
%   Examples:
%      v = loadvec('*.vc7');
%      showf(filterf(v));
%      showf(filterf(v,2,'flat'));
%      showf(filterf(vec2scal(v,'rot'),2));
%
%   See also SHOWF, MEDIANF, BWFILTERF, ADDNOISEF, INTERPF,
%   GAUSSMAT, CONV2.


%   F. Moisy, moisy_at_fast.u-psud.fr
%   Revision: 1.60,  Date: 2013/02/22
%   This function is part of the PIVMat Toolbox


% History:
% 2003/??/??: v1.00, first version.
% 2004/04/30: v1.01.
% 2004/06/26: v1.02, changed help text.
% 2005/02/23: v1.03, cosmetics.
% 2005/04/28: v1.04, propagates the field 'history'
% 2005/09/03: v1.05, arg check added.
% 2005/09/28: v1.06, option 'same' added.
% 2005/10/11: v1.10, now called filterf (from filtervec v1.06, 2005/09/28).
% 2005/10/29: v1.11, bug fixed when F is a scalar field (truncf).
% 2006/03/10: v1.12, new history
% 2006/04/28: v1.13, default value fsize=1, and accepts numeric input
% 2006/09/12: v1.20, now works with varargin; new syntax for the options.
%                    new method 'igauss' (A. Muller)
% 2008/04/15: v1.21, bug fixed in field 'history'
% 2008/09/16: v1.22, another bug fixed in field 'history'
% 2008/10/08: v1.30, missing data are now properly excluded from the
%                    convolution, using NaNs. option 'zeros' removed. use INTERPF instead.
% 2009/03/17: v1.40, FSIZE can be an array
% 2009/04/13: v1.50, main loop 'parfor' (parallel computing toolbox)
% 2010/05/03: v1.51, parfor loop removed (it is slower using parfor!)
% 2013/02/22: v1.60, works with 3D fields

error(nargchk(1,inf,nargin));

if (ischar(v) || iscellstr(v) || isnumeric(v)),
    v = loadvec(v);
end

if nargin<2
    filtsize=1;
end

if nargin<3
    method='gauss';
end

if numel(filtsize)>1
    if numel(filtsize)~=numel(v)
        error('Dimensions of vector/scalar fields F and filter sizes must be equal.');
    end
else
    filtsize = filtsize * ones(1,numel(v));
end

v=zerotonanfield(v);
fv=v;

vecmode = isfield(v(1),'vx');   % 1 for vector field, 0 for scalar field

for i=1:numel(v)    % no parfor here! (otherwise it is slower)

    if filtsize(i)~=0
        
        % defines the convolution matrix:
        if strncmpi(method,'flat',1)
            mat=ones(ceil(filtsize(i)))/ceil(filtsize(i))^2;
        elseif strncmpi(method,'gauss',2)
            mat=gaussmat(filtsize(i),1+2*ceil(3.5*filtsize(i)),'igauss');
        else
            mat=gaussmat(filtsize(i),1+2*ceil(3.5*filtsize(i)));
        end

        % computes the filtered fields:
        if any(strncmpi(varargin,'same',1))
            if vecmode
                fv(i).vx=conv2(v(i).vx,mat,'same');
                fv(i).vy=conv2(v(i).vy,mat,'same');
                if isfield(v(i),'vz')
                    fv(i).vz=conv2(v(i).vz,mat,'same');
                end
            else
                fv(i).w=conv2(v(i).w,mat,'same');
            end
        else
            if vecmode,
                fv(i).vx=conv2(v(i).vx,mat,'valid');
                fv(i).vy=conv2(v(i).vy,mat,'valid');
                if isfield(v(i),'vz')
                    fv(i).vz=conv2(v(i).vz,mat,'valid');
                end
                ntrunc=floor((size(v(i).vx,1)-size(fv(i).vx,1))/2);
            else
                fv(i).w=conv2(v(i).w,mat,'valid');
                ntrunc=floor((size(v(i).w,1)-size(fv(i).w,1))/2);
            end
            % computes the new axis (smaller than the original ones,
            %  because the borders are excluded from the convolution):
            fv(i).x = v(i).x((1+ntrunc):(length(v(i).x)-ntrunc));
            fv(i).y = v(i).y((1+ntrunc):(length(v(i).y)-ntrunc));
        end

        fv(i).history = {v(i).history{:} ['filterf(ans, ' num2str(filtsize(i)) ', ''' method ''', ' varargin{:} ')']}';

    end
end

fv=nantozerofield(fv);

if nargout==0
    showf(fv);
    clear fv
end
