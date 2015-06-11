function dout = hotwireData(PUNITS)

    %Hard-coded parameters which may need to be changed
D = 0.0254*1.5; %jet diameter - m
Pa = 100800;  %Pa - ambient pressure
	%Code assumes isentropically/ideally expanded flow
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  FILE SELECTION/BASIC CONFIGURATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Asks user for a directory containing data - can only process one
    %directory per program execution.
dir_name = uigetdir('My Computer','Specify Directory Containing Data'); 
flist = struct2cell(dir(dir_name));
q = strfind(flist(1,:),'.DATA'); %extracts list of data files - ignores all other files
keep = zeros(size(q));
for n = 1:length(keep)
    if ~isempty(q{n})
        keep(n) = 1;
    end
end
flist = flist(1,logical(keep));
[keep,ok] = listdlg('PromptString','Select files to process:',...
                'SelectionMode','multiple',...
                'ListString',flist);    %Asks the user to select the subset of data files they wish to process
if ok==0
    error('PT_UserCancel','Program Terminated Due to User Selection of Cancel')
end
sampleRate = input('  Enter the Sample Rate for Data: ');   %Asks the user for the rate at which the data was acquired - If user enters nothing (Press Enter), program defaults to 200 kSamples/sec
if isempty(sampleRate)
    sampleRate = 200000;
    disp('    Sample Rate Set to Default: 200000')
end
dSize = input('  Enter the number of data points in a block: ');    %Asks the user for the number of data points in a block on which FFT will be calculated - If user enters nothing (Press Enter), program defaults to 8192 points
if isempty(dSize)
    dSize = 8192;
    disp('    Block Size Set to Default: 8192');
end

DPoints = 1;
fx = (0:ceil((dSize+1)/2)-1)*sampleRate/dSize;  %Calculates frequency axis
fx = fx(DPoints:end-DPoints); %First and last points will be thrown out as garbage

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  DATA PROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %All operations beyond this point are done without user interaction.
    
POS = [];   %Location of measurement - um
MV = [];    %mean voltage - V
avgFFT{1} = [];

M = [];
Temp = [];
U = [];
ReN = [];
for n = 1:length(keep)
    Fnm = flist{keep(n)};
    disp(['    Processing File: ' Fnm])  %Displays file currently being processed
        %Locates pressure information in file name
    PLoc = 0;   
    Pend = min(strfind(Fnm(PLoc+1:end),'_'))+PLoc;
    if isempty(Pend)
        Pend = length(Fnm)-length('.DATA')+1;
    end
    Ptxt = Fnm(PLoc+2:Pend-1);
    Ptxt(strfind(Ptxt,'p'))='.';    %if p was used instead of . in floating point number, it is replaced
    Po(n) = str2num(Ptxt);
    if strmatch(PUNITS,'inWater')
        Po(n) = Po(n)*248.9;    %Convert to Pa gauge pressure
    elseif strmatch(PUNITS,'PSIG')
        Po(n) = Po(n)*6.8927e3;
    else
        error('No Units')
    end
        %Locates temperature information in file name
    TLoc = strfind(Fnm,'_T');   
    Tend = min(strfind(Fnm(TLoc+1:end),'_'))+TLoc;
    if isempty(Tend)
        Tend = length(Fnm)-length('.DATA')+1;
    end
    Ttxt = Fnm(TLoc+2:Tend-1);
    Ttxt(strfind(Ttxt,'p'))='.';    %if p was used instead of . in floating point number, it is replaced
    To = str2num(Ttxt)+273.15;  %Converts to Kelvin
        %Locates Probe position information in file name
    FLoc = strfind(Fnm,'_Pos');   
    if ~isempty(FLoc)
        Fend = min(strfind(Fnm(FLoc+1:end),'_'))+FLoc;
        if isempty(Fend)
            Fend = length(Fnm)-length('.DATA')+1;
        end
        POS(n) = str2num(Fnm(FLoc+4:Fend-1));
    else
        error('No Position');
    end
    
    M(n) = sqrt((((Po(n)+Pa)/Pa)^((1.4-1)/1.4) -1)*2/(1.4-1));    %Mach Number
    Temp(n) = To/(1+1/5*M(n)^2);  %Calculates isentropically expanded jet exit temperature - K
    U(n) = M(n)*sqrt(7/5*287.05*Temp(n));   %velocity - m/s
    ReN(n) = 1.01e5/(287.05*Temp(n)) *U(n)*...
        D /(1.827e-5*(0.555*291.15+66.667)/(0.555*Temp(n)+66.667).*(Temp(n)/291.15).^(3/2));  %Calculates jet exit Reynolds number using Sutherlands formula for viscosity and ideal gas law for density
    
    
    rRaw = dlmread([dir_name filesep flist{keep(n)}],'\t'); %reads data file
    MV(n) = mean(rRaw); %average voltage - V
    
    S = size(rRaw);
    nBlocks = S(1)/dSize;   %Determines number of blocks in file - Note: All channels are processed simultaneously
    rData = rRaw(1:dSize,:)-repmat(mean(rRaw(1:dSize,:),1),[dSize,1]);  %Subtracts mean value from first block for every column
    rFFT = fft(rData,dSize,1);  %calculates FFT of first block for all channels 
    for m = 2:nBlocks
        rData = rRaw(dSize*(m-1)+1:dSize*m,:)-repmat(mean(rRaw(dSize*(m-1)+1:dSize*m,:),1),[dSize,1]);
        rFFT(:,:,m) = fft(rData,dSize,1);
    end
    rFFT = rFFT(1:ceil((dSize+1)/2),:,:);   %Discards symmetric portion
    mFFT = rFFT.*conj(rFFT)/dSize;  %Scales spectrum according to block size 
    if rem(dSize,2) %Conserves energy which would otherwise be lost by discardiing symmetric portion
        mFFT(2:end,:,:) = mFFT(2:end,:,:)*4;
    else
        mFFT(2:end-1,:,:) = mFFT(2:end-1,:,:)*4;
    end
    avgFFT{n} = sum(mFFT,3)/nBlocks/U(n)^2; %Calculates MS average spectrum and scales by local velocity
end

dout.Po = Po;   %Stagnation pressure (psig)

[dout.POS,IX] = sort(POS);
dout.M = M(IX);
dout.Temp = Temp(IX);
dout.MV = MV(IX);
dout.avgFFT = avgFFT(IX);
dout.fx = fx;
dout.U = U(IX);
dout.ReN = ReN(IX);

