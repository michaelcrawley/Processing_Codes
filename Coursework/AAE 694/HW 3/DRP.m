function [coefs] = DRP(N,h,varargin)
%Calculate coefficients for use in finite difference scheme using
%Dispersion-Relation-Preserving schemes (central differencing only)
%following code only works for first order derivatives
%Completed by Michael Crawley for AAE694 HW3

    D = 1; %Code currently only works for first order derivatives
    A = zeros(2*N+1,2*N+1); %initialize left matrix for solution (Ax=F)
    F = zeros(2*N+1,1); %initialize right matrix for solution
    F(D+1) = factorial(D); %Set derivative dependence
    if isempty(varargin)
        PPW = 4;
    else
        PPW = varargin{1};
    end
    delta = 2*pi/PPW;
    
    %Set N:N-2 rows in left matrix based off of taylor series expansions
    for j=1:2*N-1
        A(j,:) = (-N:N).^(j-1);
    end

    %Set second to last row based off of left-right taylor series property
    A(2*N,:) = [(-1)^(D+1)*ones(1,N) 0 ones(1,N)];

    %calculate d/dan coefficients
    da = [1;A(1:2*N,2:2*N+1)\-A(1:2*N,1)];

    %Create integration string
    str = ['(',num2str(da(1)),'*i*exp(i*',num2str(-N),'*b)'];
    for j = -N+1:N
        str = strcat(str,['+(',num2str(da(j+N+1)),')*i*exp(i*',num2str(j),'*b)']);
    end
    str = [str,')'];

    %Integrate for final linear equation
    Y = double(int(['-b*',str],-delta,delta)); %Why does this need to be -b instead of b?
    F(end) = Y;
    for j = -N:N
        A(N+N+1,j+N+1)= double(int(['i*exp(i*',num2str(j),'*b)*',str],-delta,delta));
    end

    %Calculate coefficients from Ax=F
    coefs = (A\F)/(h^D);
end

