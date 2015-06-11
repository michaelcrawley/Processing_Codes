%The purpose of this file is to serve as a place holder for temporary
%matlab scripts

NFCh = 12:19;
Lch = length(NFCh);
D = 0.0254;
dx = 0.5*D;
dy = 0.2*D;
dt = 1/(2e5);
a = sqrt(1.4*287*(25.6+273.1));
theta = 8; %array angle
transform = [secd(theta) -tand(theta); 0 1]; %inversion matrix
windowfun = @(x) tukeywin(x,0.1);
funname = 'mod-tukey cabana';

for n = 1:length(files)
    load([src filesep files{n}]);
    [~,fname] = fileparts(files{n});
    Lsm = min(cellfun(@(x) length(x), pf.sm_wvfm)); %Get minimum length for averaged, filtered waveform
    Lf = length(pf.sm_wvfm);

    %Pull waveform Data
    waveform = zeros(Lsm,Lch,Lf);
    X = zeros(Lf,Lch);
    Y = X;
    for nn = 1:Lf
       waveform(:,:,nn) = pf.sm_wvfm{nn}(1:Lsm,NFCh);
       X(nn,:) = pf.phys{nn}.x;
       Y(nn,:) = pf.phys{nn}.y;
    end
    waveform = permute(waveform,[3 2 1]);

    %Reshape into correct grid
    [X Y waveform] = ReshapeGrid(X,Y,waveform);
    [M,N,P] = size(waveform);
    xm = mean(mean(mean(waveform)));

    %Perform decomposition
    [hydro acoustic f PSD S k vel subW supW] = Phase_Velocity_Decompose(waveform,[dy dx dt],a,transform,windowfun);
      
    %Compute Prms and OASPL maps
    pref = 2e-5;
    hmean = repmat(mean(hydro,3),[1 1 P]);
    prms.hydro = sqrt(mean((hydro-hmean).^2,3));
    prms.acoustic = sqrt(mean(acoustic.^2,3));
    prms.waveform = sqrt(mean(waveform.^2,3));
    OASPL.hydro = 20*log10(prms.hydro/pref);
    OASPL.acoustic = 20*log10(prms.acoustic/pref);
    OASPL.waveform = 20*log10(prms.waveform/pref);

    %Save data
    save([pwd filesep fname ' decomp ' funname '.mat'],'waveform','X','Y','PSD','f','S','xm','k','vel','subW','supW','hydro','acoustic','prms','OASPL','windowfun');
end