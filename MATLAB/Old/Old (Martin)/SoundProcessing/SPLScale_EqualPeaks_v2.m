%This function reads in the spectrum data created in SoundProcess.m and
%performs a SPL correction.  The correction is force all peaks to coincide with a value of zero.
%Output is a series of 
%.fig and .png files displaying the scaled spectra in a directory chosen by 
%the user. Files created are:
%	Ch(N)_SmoothCompare - displays channel N spectra after moving average has been applied to the spectra (smoothed)
%	Ch(N)_Scale - displays channel N spectra scaled to the reference - values of n for each spectrum are printed on the top of the graph
%	Ch(N)_SmoothScale - displays channel N scaled, smoothed spectra

clear variables
clc
disp('All Spectra Must Have Same Number of Channels in Same Order')
numSp=input('Enter the number of spectra to be compared: ');


done2 = false;
toto = 0;
n = 0;
SPfile = [];
M = [];
M2 = [];
TempRatio = zeros(numSp,1);
JetVelRatio = zeros(numSp,1);
Re = zeros(numSp,1);
Mag = zeros(numSp,2);
while ~done2
    SPpath = uigetdir('My Computer','Specify Directory Containing Data');
    flist = struct2cell(dir(SPpath));
    q = strfind(flist(1,:),'.fftNOS');
    keep = zeros(size(q));
    for n = 1:length(keep)
        if ~isempty(q{n})
            keep(n) = 1;
        end
    end
    flist = flist(1,logical(keep));
    [keep,ok] = listdlg('PromptString','Select files to process:',...
                    'SelectionMode','multiple',...
                    'ListString',flist);
    if ok==0
        error('Program Terminated Due to User Selection of Cancel')
    end
    SPf = flist(keep);
    

    SPfile = [SPfile SPf];
    totn = toto+length(SPf);
    for n = toto+1:totn
        M=dlmread([SPpath filesep SPfile{n}],'\t',1,0);
        if M(1,1) == 0
            M = M(2:end,:);
        end
        if exist([SPpath filesep SPfile{n}(1:end-6) 'SET'],'file')
            done = false;
            fid = fopen([SPpath filesep SPfile{n}(1:end-6) 'SET'],'r');
            while ~done
                L = fgetl(fid);
                mtch = ~isempty(strmatch('Jet Temperature Ratio:',L));
                if mtch
                    TempRatio(n) = str2num(L(24:end));
                    done = true;
                end
            end
            done = false;
            while ~done
                L = fgetl(fid);
                mtch = ~isempty(strmatch('Jet Velocity Ratio:',L));
                if mtch
                    JetVelRatio(n) = str2num(L(21:end));
                    done = true;
                end
            end
            done = false;
            while ~done
                L = fgetl(fid);
                mtch = ~isempty(strmatch('Reynolds Number:',L));
                if mtch
                    Re(n) = str2num(L(18:end));
                    done = true;
                end
            end
            done = false;
            while ~done
                L = fgetl(fid);
                mtch = ~isempty(strmatch('  Magnitude (dB):',L));
                if mtch
                    Mag(n,:) = str2num(L(18:end));
                    done = true;
                end
            end      
            fclose(fid);
        elseif exist([SPpath filesep SPfile{n}(1:end-9) 'SET'],'file')
            done = false;
            fid = fopen([SPpath filesep SPfile{n}(1:end-9) 'SET'],'r');
            while ~done
                L = fgetl(fid);
                mtch = ~isempty(strmatch('Jet Temperature Ratio:',L));
                if mtch
                    TempRatio(n) = str2num(L(24:end));
                    done = true;
                end
            end
            done = false;
            while ~done
                L = fgetl(fid);
                mtch = ~isempty(strmatch('Jet Velocity Ratio:',L));
                if mtch
                    JetVelRatio(n) = str2num(L(21:end));
                    done = true;
                end
            end
            done = false;
            while ~done
                L = fgetl(fid);
                mtch = ~isempty(strmatch('Reynolds Number:',L));
                if mtch
                    Re(n) = str2num(L(18:end));
                    done = true;
                end
            end
            done = false;
            while ~done
                L = fgetl(fid);
                mtch = ~isempty(strmatch('  Magnitude (dB):',L));
                if mtch
                    Mag(n,:) = str2num(L(18:end));
                    done = true;
                end
            end
            fclose(fid);
        else
            error('Dummy');
%             disp(['  File ' num2str(1) ': ' SPfile ' Info - '])
%             TempRatio(1) = input('Input temperature ratio (found in the SET file) for the spectrum: ');
%             JetVelRatio(1) = input('Input jet velocity ratio (found in the SET file) for the spectrum: ');
        end
        S = size(M);
        Nch = S(2)-2;
        M(:,3:Nch+2) = M(:,3:Nch+2) - repmat(Mag(n,:),length(M(:,1)),1);
        M2{n} = M;
               
    end
    if totn >= numSp
        done2 = true;
        numSp = totn;
    else
        toto = totn;
        disp(['Select Next Spectra: ' num2str(toto+1)])
    end
end

disp('Select Output Directory')
SvPath = uigetdir('My Computer','Select Output Directory');

colors = 'krbgmc';
markers = '.ox+*sdv^<>ph';
for n = 1:Nch
    %Plot Scaled Ch(s)
    clear xx yy mrk;
    figure
    specNameM=[];
    for m=1:numSp
        xx{m} = M2{m}(:,2);
        yy{m} = M2{m}(:,n+Nch);
        mrk{m,1} = markers(mod(m,13)+1);
        mrk{m,2} = colors(mod(m,6)+1);
        specName{m}=num2str(TempRatio(m),'%.2f');
    end
    semilogxSparseMarker(xx,yy,mrk,5,specName)
    grid on
    xlabel('St_D')
    ylabel('SPL (dB)')
    title(['Scaled Channel ' num2str(n) ' SPL Data'])
    saveas(gcf,[SvPath filesep 'Ch' num2str(n) '_Scaled.png'])
    saveas(gcf,[SvPath filesep 'Ch' num2str(n) '_Scaled.fig'])
    close
    
end
save([SvPath filesep 'Composite.mat'],'JetVelRatio','M2','Mag','TempRatio','Re')
    

