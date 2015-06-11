function [rec] = iSpectralPODnD(alpha,phi)
%This function performs a reconstruction of POD modes and coefficients. The
%alpha and phi matrices need to be organized per the output of
%SpectralPODnD. By default, this function will perform the reconstruction
%over all of the POD modes (second index in alpha and phi). In order to
%perform a low-dimensional reconstruction, the unwanted modes must either
%be set to zero or removed from the matrices entirely.

    [NF,NM,NB] = size(alpha);
    NX = size(phi,1);
    rec = complex(zeros(NF,NX,NB));
    for n1 = 1:NF %frequency - looped
        for n3 = 1:NB %Block number - looped
            for n2 = 1:NM %modes will be summed
                rec(n1,:,n3) = rec(n1,:,n3) + alpha(n1,n2,n3)*phi(:,n2,n1)';
            end
        end
    end
    
    rec = real(ifft(rec,[],1));
end