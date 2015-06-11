function [recon,A] = TemporalLSE(field,signal) 
%This function computes a least-squares correlation matrix between a given
%field and a reference signal, per the theory outlined by Adrian regarding
%conditional averaging based on linear stochastic estimation. This code
%computes the estimate using single-time, zero-time lag cross-correlations.
%The field and reference signals need to be ordered as 2-D
%matrices, the the first dimension being time and the second being
%space/component/etc.
%Inputs:
%           field:      Generalized field for conditional averaging
%           signal:     Generalized reference signal
%Outputs:
%           recon:      Reconstructed field per signal
%           A:          Correlation matrix
%Last updated: 2015-03-04 by Michael Crawley

    unconditional = signal.'*signal;
    conditional = signal.'*field;
    A = unconditional\conditional;
    recon = signal*A;
end