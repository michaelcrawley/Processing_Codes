function WriteLAFPAWaveFormFile(F,D,M,M2,N,ModF,ModA,rotFlag,RA,CH,RForDC)

    clockpulse = 4000000;
    
    %check inputs
    if ~exist('rotFlag','var'), rotFlag = 0; end
    if ~exist('RA','var'), RA = 0; end
    if ~exist('CH','var'), CH = ones(N); end
    if ~exist('RForDC','var'), RForDC = 0; end
    
    [PS, NP] = calcPS(F,D,M,M2,N,rotFlag,RA,CH,clockpulse);
    [PSt,NPt] = Modulate(ModF,ModA,NP,PS,clockpulse);
    [PA,~] = calcPA(PSt, N, NPt,CH,RForDC);
    
    filename = ['F',num2str(F),'_M',num2str(M),'_M',num2str(M2),'_PW',num2str(D),'_ModF',num2str(ModF),'_ModA',num2str(ModA),'.lwv'];
    
    fid = fopen(filename,'w');
    fprintf(fid,'%d\t',PA);
    fclose(fid);


end