%The purpose of this file is to serve as a place holder for temporary
%matlab scripts
%This function performs a 2-D spatial Fourier transform on the
%phase-averaged nearfield pressure waveforms

NFCh = 12:19;
FFCh = 1:12;
Lch = length(NFCh);
D = 0.0254;
dx = 0.5*D;
dy = 0.2*D;
dt = 1/(2e5);
a = sqrt(1.4*287*(25.6+273.1));
theta = 8; %array angle
windowfun = @(x) tukeywin(x,0.1);

for n = 1:length(files)
    load([src filesep files{n}]);
    [~,fname] = fileparts(files{n});
    Lsm = min(cellfun(@(x) length(x), pf.sm_wvfm)); %Get minimum length for averaged, filtered waveform
    Lf = length(pf.sm_wvfm);

    %Pull waveform Data
    waveform = zeros(Lsm,Lch,Lf);
    FFwaveform = zeros(Lsm,length(FFCh));
    X = zeros(Lf,Lch);
    Y = X;
    for nn = 1:Lf
       waveform(:,:,nn) = pf.sm_wvfm{nn}(1:Lsm,NFCh);
       FFwaveform = FFwaveform + pf.sm_wvfm{nn}(1:Lsm,FFCh)/Lf;
       X(nn,:) = pf.phys{nn}.x;
       Y(nn,:) = pf.phys{nn}.y;
    end
    waveform = permute(waveform,[3 2 1]);

    %Reshape into correct grid
    [X Y waveform] = ReshapeGrid(X,Y,waveform);
    [M,N,P] = size(waveform);
    xm = mean(mean(mean(waveform)));
    
    %Perform Spatial Fourier Transform
    [PSD f] = PSDN(waveform,[1 2],[dy dx],windowfun);
    PSD = fftshift(PSD); 
        
    %Change to orthogonal angular wavenumbers
    angf = cellfun(@(x) fftshift(2*pi*x),f,'uniformoutput',false);
    angf{2} = angf{2}*cosd(theta); %transform from km to kx
    
    %Save data
    save([pwd filesep fname ' 2DspatialFFT.mat'],'waveform','FFwaveform','X','Y','angf','PSD','windowfun');
end