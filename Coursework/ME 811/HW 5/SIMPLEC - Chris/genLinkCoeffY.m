function [B] = genLinkCoeffY(BCu,BCv,u,v,p,rho,mu)
% Calculate the link coefficients for the y-momentum equation.

%% Transpose the inputs
BC1.t = BCv.r';		BC1.r = BCv.t';
BC1.b = BCv.l';		BC1.l = BCv.b';

BC2.t = BCu.r';		BC2.r = BCu.t';
BC2.b = BCu.l';		BC2.l = BCu.b';

%% Call the already formulated x-momentum link coefficient generator
[A] = genLinkCoeffX(BC1,BC2,v',u',p',rho,mu);

%% Transpose the output
B.s = A.w';			B.w = A.s';
B.e = A.n';			B.n = A.e';
B.o = A.o';			B.p = A.p';