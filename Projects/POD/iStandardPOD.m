function [rec] = iStandardPOD(phi,ak)
%Phi is the eigenvectors (each column is a different POD mode)
%ak are the time-dependent coefficients
    rec = phi*ak;

end