function [R phase f] = computeXSpectra(src_dir,flist)
%this function computes the cross-spectral density matrix 


    %Temporary constants
    FS = 2e5;
    
    %load data file 
    load([src_dir filesep flist]);
    
    %Get size info
    [BS,Ch,~] = size(nf.pblocks.p);
    
    %initialize variables
    R = zeros(ceil((BS+1)/2),Ch,Ch);
    phase = R;
        
    for n = 1:Ch
        for nn = 1:Ch
            [~,f,phase(:,n,nn),r] = coherence(squeeze(nf.pblocks.p(:,n,:)),squeeze(nf.pblocks.p(:,nn,:)),FS,BS,@hann);
            R(:,n,nn) = real(r);
        end        
    end


end