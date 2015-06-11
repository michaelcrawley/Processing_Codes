function [raw,subsonic,supersonic,f] = computePSD_decomposed(src_dir,flist)
%This only works on single files

    tmp = load([src_dir filesep flist{1}]);
    
    %Reshape
    N = size(tmp.data.intwaveform);
    raw = reshape(permute(tmp.data.intwaveform,[1 3 2]),[N(1)*N(3),N(2)]);
    subsonic = reshape(permute(tmp.data.subsonic,[1 3 2]),[N(1)*N(3),N(2)]);
    supersonic = reshape(permute(tmp.data.supersonic,[1 3 2]),[N(1)*N(3),N(2)]);
    
    %Compute PSD
    [PSD,f] = calSPL([raw,subsonic,supersonic],N(1),tmp.daq_params.phys.FS,'hamming',0,[0 0]);
    
    %Extract outputs
    raw = PSD(:,1:N(2));
    subsonic = PSD(:,(N(2)+1):2*N(2));
    supersonic = PSD(:,(2*N(2)+1):end);
end