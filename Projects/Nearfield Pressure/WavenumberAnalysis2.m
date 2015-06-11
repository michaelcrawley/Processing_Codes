%This code performs a 2-D Fourier transform in the x,t domain on the
%filtered, non-phase-averaged nearfield pressure data

%Constants
D = 0.0254;
dx = D;
dt = 1/(2e5);
theta = 8;
windowfun = @(x) tukeywin(x,0.1);

for n = 1:length(files)
    %load data
    load([src filesep files{n}]);
    [~,fname] = fileparts(files{n});
    
    %Get Data info
    [BS,Lch,NB] = size(nf.pblocks.smp);
    
    %Calc kx,w PSD
    [PSD f] = PSDN(nf.pblocks.smp,[1 2],[dt dx],windowfun);        
    avgPSD = mean(PSD,3);
    avgPSD = fftshift(avgPSD);
    
    %Convert to angular and from km to kx
    angf = cellfun(@(x) fftshift(2*pi*x),f,'uniformoutput',false);
    angf{2} = angf{2}*cosd(theta); %transform from km to kx
    
    %Save data
    save([save_dir filesep fname '.mat'],'src_dir','phys','windowfun','filename','PSD','avgPSD','f','angf');  
end