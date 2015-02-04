function V = solveY(AO,AE,AW,AN,AS,S,Vold,sweeps)

%this solves the x-momentum equation when given the coefficients, a guess
%velocity field, and the number of sweeps of ADI to be performed

[N,M] = size(AO);
K = M*N;

AO = reshape(AO',K,1);
AE = reshape(AE',K,1);
AW = reshape(AW',K,1);
AN = reshape(AN',K,1);
AS = reshape(AS',K,1);
S = reshape(S',K,1);

[V,~] = ADI(AS,AW,AO,AE,AN,S,Vold,num2str(sweeps));