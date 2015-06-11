function [F,OA,NTOA,TE] = ToneEnergyContribution(SD)

Res = 200000/8192;
SD = 'F:\Research\Mach0.9\ForcingNoiseData\Array Measurements\Recollecting 200909';

    %Asks user for a directory containing data - can only process one
    %directory per program execution.
dir_name = uigetdir(SD,'Specify Directory Containing Data'); 
flist = struct2cell(dir(dir_name));
q = strfind(flist(1,:),'.fftNOS'); %extracts list of data files - ignores all other files
keep = ~cellfun('isempty',q);
flist = flist(1,logical(keep));
q = strfind(flist(1,:),'.S.fftNOS'); %extracts list of data files - ignores all other files
keep = cellfun('isempty',q);
flist = flist(logical(keep));
[keep,ok] = listdlg('PromptString','Select files to process:',...
                'SelectionMode','multiple',...
                'ListString',flist);    %Asks the user to select the subset of data files they wish to process
if ok==0
    error('Program Terminated Due to User Selection of Cancel')
end
flist = flist(keep);

for n = 1:length(flist)
    Fnm = flist{n};
    FLoc = strfind(Fnm,'_F');   
    if ~isempty(FLoc)
        Fend = min(strfind(Fnm(FLoc+1:end),'_'))+FLoc;
        F(n) = str2num(Fnm(FLoc+2:Fend-1))*1000; clear FLoc Fend;
    else
        error('Not a forced case')
    end

    d = dlmread(fullfile(dir_name,Fnm),'\t',1,0);
    Bg = min(find(d(:,2) > 0.04)); %Excludes spectral components below 0.04 Strouhal number
    Ed = min(find(d(:,2) > 4));    %Excludes spectral components above 4 Strouhal number
    x = d(Bg:Ed,1); d = d(Bg:Ed,3:end);
    
    fm = (1:1:max(x)/F(n))*F(n);  %Forcing harmonic frequencies
    fi = round((fm-x(1)+Res)/Res);  %calculates index of harmonics based on sampling parameters
    fi = fi(fi < Ed-Bg-4+1);   %Excludes points beyond integration range
    
    OA(n,:) = sum(10.^(d/10),1)*Res; %Calculates OASPL of signal with tones - (square pressure ratio units)
    il = true(size(x));   %initializes logical for integration
    IL = zeros(size(fi)); IM = IL; REP = zeros(1,size(d,2)); T = zeros(1,size(d,2));    %Initialization
    for m = 1:length(fi)
        [C,I] = max(d(fi(m)-4:fi(m)+4,:),[],1);  %Finds peak of harmonic across all channels
        I = mode(I);    %Uses most common index as reference for all channels
        I = I + fi(m)-4 -1; IL(m) = I-2;   IM(m) = I+2;   %calculates harmonic point range [IL,IM]
        il(IL(m):IM(m)) = false;
        T = T + sum(10.^(d(IL(m):IM(m),:)/10),1)*Res; %accumulates the energy contained the harmonics
        REP = REP + 10.^((d(IL(m),:)+d(IM(m),:))/20)*Res*(IM(m)-IL(m)+1);  %calculates tone-removed replacement level
    end
    TE(n,:) = T;
    NTOA(n,:) = sum(10.^(d(il,:)/10),1)*Res +REP;
%     semilogx(x(il),d(il,1),'LineWidth',2)
%     hold on
%     semilogx(x,d(:,1),'r')
%     title(num2str(F(n)))
%     hold off
end
    
[F,IX] = sort(F);
OA = OA(IX,:);
NTOA = NTOA(IX,:);
TE = TE(IX,:);
