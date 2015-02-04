function [A] = genLinkCoeffT(bc,u,v,p,rhof,cpf,kf,M,N,m,n,dy,dx)

rcf.x = rhof.x .* cpf.x;
rcf.y = rhof.y .* cpf.y;

A = genLinkCoeffMom(bc,u,v,p,rcf,kf,M,N,m,n,dy,dx);

%% Calculate the source links
A.p = zeros(M,N);		% Overwrites A.p.x and A.p.y from genLinkCoeffMom()

% Step top
i=m+1; j=1:n;
A.p(i,j) = kf.y(i,j).*bc.stop.T;

% Step face
i=1:m; j=n+1;
A.p(i,j) = kf.x(i,j).*bc.sface.T;