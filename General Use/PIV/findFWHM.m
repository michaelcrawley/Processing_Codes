function fwhm = findFWHM(x,y,U,varargin)
%Finds the full width at half maximum of mean velocity using the provided
%velocity component. If varargin in isn't empty, program assumes the
%additional argument is the width cutoff parameter (i.e. width at 0.7 of
%max). Code assumes y-coordinates are descending.

if isempty(varargin)
    L = 0.5;
else
    L = varargin{1};
end

fwhm = zeros(size(x(:,1)));
for n = 1:length(x(:,1))    
    Uf = U(n,:);
    [C,I] = max(Uf);
    Uf = Uf/C;  %Normalize the flow profile by the local peak velocity
    
    if C > 0    %Only look for percent-max if the peak has positive value
        yp = max(find(Uf(1:I) < L));    %Find percent-max to left of centerline
        if isempty(yp)  %If percent-max condition not satisfied, use left end of domain
            yp = 1;     %yp==1 is a trigger for a later event
            ypx = y(1,1);
        else    %Otherwise, interpolate points near percent-max to find exact location
            ypx = linspace(y(1,yp),y(1,yp+1),100);  
            Ufi = interp1(y(1,yp:yp+1),Uf(yp:yp+1),ypx);
            yp = max(find(Ufi < L));
        end

        ym = min(find(Uf(I:end) < L))+I-1;  %Find percent-max to right of centerline
        if isempty(ym)  %If percent-max condition not satisfied, use right end of domain
            ym = 1;     %ym==1 is a trigger for a later event
            ymx = y(1,end);
        else    %Otherwise, interpolate points near percent-max to find exact location
            ymx = linspace(y(1,ym-1),y(1,ym),100);
            Ufi = interp1(y(1,ym-1:ym),Uf(ym-1:ym),ymx);
            ym = min(find(Ufi < L));
        end
        if yp==1 & ym~=1    %If left side of domain was not determined, double the right side value
            fwhm(n) = -2*ymx(ym);
        elseif yp~=1 & ym==1    %If right side of domain was not determined, double the left side value
            fwhm(n) = 2*ypx(yp);
        else    %If both sides were missing, assume width is the full extent of the domain
            fwhm(n) = ypx(yp)-ymx(ym);
        end
    else
        fwhm(n) = 0;    %In the event of no-local max, assign zero.
    end
end

    %Replaces zeros (assumed to be at the ends of the profile) with the
    %nearest non-zero neighbor.
Q = find(fwhm==0); Q = Q(:)';
cp = round(length(x(:,1))/2);
Ql = fliplr(Q(Q < cp));
for n = 1:length(Ql)
    fwhm(Ql(n)) = fwhm(Ql(n)+1);
end
Ql = Q(Q >= cp);
for n = 1:length(Ql)
    fwhm(Ql(n)) = fwhm(Ql(n)-1);
end

    %A high order polynomial fit to the data smooths out any
    %inconsistencies created by calculating width on half or whole domain.
warning('off','MATLAB:polyfit:RepeatedPointsOrRescale');
Cp = polyfit(x(:,1),fwhm,10);
ft = polyval(Cp,x(:,1));

d = fwhm-ft;
rp = find(d > mean(d)+3*std(d));
rp = rp(rp>1); rp = rp(rp<length(d));
for n = 1:length(rp)
    fwhm(rp(n)) = (fwhm(rp(n)+1)+fwhm(rp(n)-1))/2;
end
