function [pp] = ProcessBaselineResults(pp,out_dir)
%Finds processed baseline .SET data files, and reads in the results for use
%in calculating dOASPL, dAAE, etc for the forced cases.  This function must
%be run after all the baseline files are processed but before any of the
%forced cases are processed.  Baseline results are added to the processing
%parameters structure.

%Last updated by Michael Crawley on 2011-07-04

    slist = struct2cell(dir(out_dir));
    qs = strfind(slist(1,:),'Baseline');
    qs2 = strfind(slist(1,:),'.SET');
    kp = zeros(size(qs));
    for m = 1:length(qs)
        if (~isempty(qs{m}))&&(~isempty(qs2{m}))
            kp(m) = 1;
        end
    end
    if sum(kp) > 0  %If baseline SET file(s) are found
        slist = slist(1,logical(kp));
        bdata.OASPLB = zeros(1,pp.Nch);  bdata.AAEB = zeros(1,length(slist));
        bdata.OASPLB_DT = bdata.OASPLB; bdata.AAEB_DT = bdata.AAEB; bdata.BT = bdata.AAEB; bdata.NJPB = bdata.AAEB; bdata.NJPB_DT = bdata.AAEB;

        if ~exist('SETpnts','var')
                %Gets the SET file line numbers which contain the needed
                %information.  
            q = strfind(slist{1},'_');
            if ~isempty(q)
                fp{1} = slist{1}(1:q(1)-1);
                if length(q) > 1
                    for nn = 2:length(q)
                        fp{nn} = slist{1}(q(nn-1)+1:q(nn)-1);
                    end
                else
                    nn = 1;
                end
                fp{nn+1} = slist{1}(q(end)+1:end);

                sOffset = 0;
                for nn = 1:length(fp)
                    sOffset = sOffset+1;
                    if strmatch(fp{nn}(1),'M')   %Mach number is already in the SET file
                        sOffset = sOffset-1;
                    elseif strmatch(fp{nn}(1),'F')   %Forcing Frequency is already in the SET file
                        sOffset = sOffset-1;
                    elseif strmatch(fp{nn}(1),'T')   %Stagnation temperature is already in the SET file
                        sOffset = sOffset-1;
                    elseif strmatch(fp{nn}(1:2),'Ba')  %Ignores Baseline tag
                        sOffset = sOffset-1;
                    end
                end
            end

            done = false; SETpnts = zeros(1,7);
            fid2 = fopen([out_dir '\' slist{1}],'r');
            while ~done
                L = fgetl(fid2);
                SETpnts(1) = SETpnts(1)+1;
                done = ~isempty(strmatch('Stagnation Temperature (K):',L));
            end
            done = false; fseek(fid2,0,'bof');
            while ~done
                L = fgetl(fid2);
                SETpnts(2) = SETpnts(2)+1;
                done = ~isempty(strmatch('OASPL (dB):',L));
            end
            done = false; fseek(fid2,0,'bof');
            while ~done
                L = fgetl(fid2);
                SETpnts(3) = SETpnts(3)+1;
                done = ~isempty(strmatch('Average Energy (dB):',L));
            end
            done = false; fseek(fid2,0,'bof');
            while ~done
                L = fgetl(fid2);
                SETpnts(4) = SETpnts(4)+1;
                done = ~isempty(strmatch('Normalized Jet Power (dB):',L));
            end
            done = false; fseek(fid2,0,'bof');
            while ~done
                L = fgetl(fid2);
                SETpnts(5) = SETpnts(5)+1;
                done = ~isempty(strmatch('OASPL_Smoothed (dB):',L));
            end
            done = false; fseek(fid2,0,'bof');
            while ~done
                L = fgetl(fid2);
                SETpnts(6) = SETpnts(6)+1;
                done = ~isempty(strmatch('Average Energy_Smoothed (dB):',L));
            end
            done = false; fseek(fid2,0,'bof');
            while ~done
                L = fgetl(fid2);
                SETpnts(7) = SETpnts(7)+1;
                done = ~isempty(strmatch('Normalized Jet Power_Smoothed (dB):',L));
            end
            fclose(fid2);
            SETpnts = SETpnts +sOffset;
        end

        for nn = 1:length(slist)    %reads baseline SET file(s) and extracts OASPL and total measured energy
            [NOASPLB,VB,SETpnts] = readSETfile([out_dir '\' slist{nn}],SETpnts);
            bdata.BT(nn) = VB{1}; %Baseline temperature (K)
            bdata.OASPLB(nn,:) = VB{2};
            bdata.AAEB(nn) = VB{3};
            bdata.NJPB(nn) = VB{4};
            bdata.OASPLB_DT(nn,:) = VB{5};
            bdata.AAEB_DT(nn) = VB{6};
            bdata.NJPB_DT(nn) = VB{7};
        end
        if length(slist) > 1
            [bdata.BT,IX] = sort(bdata.BT); %Ensures baseline data is in increasing order of temperature
            bdata.OASPLB = bdata.OASPLB(IX,:);
            bdata.OASPLB_DT = bdata.OASPLB_DT(IX,:);
            bdata.AAEB = bdata.AAEB(IX);
            bdata.AAEB_DT = bdata.AAEB_DT(IX);
            bdata.NJPB = bdata.NJPB(IX);
            bdata.NJPB_DT = bdata.NJPB_DT(IX);

            uBT = unique(bdata.BT);   %Averages data non-distinct temperatures
            if length(uBT)~=length(bdata.BT)
                Tmp = cell(6,1);
                for nn = 1:length(uBT)
                    q = find(bdata.BT==uBT(nn));
                    Tmp{1}(nn,:) = mean(bdata.OASPLB(q,:),1);
                    Tmp{2}(nn,:) = mean(bdata.OASPLB_DT(q,:),1);
                    Tmp{3}(nn) = mean(bdata.AAEB(q));
                    Tmp{4}(nn) = mean(bdata.AAEB_DT(q));
                    Tmp{5}(nn) = mean(bdata.NJPB(q));
                    Tmp{6}(nn) = mean(bdata.NJPB_DT(q));
                end
                bdata.BT = uBT;
                bdata.OASPLB = Tmp{1};
                bdata.OASPLB_DT = Tmp{2};
                bdata.AAEB = Tmp{3};
                bdata.AAEB_DT = Tmp{4};
                bdata.NJPB = Tmp{5};
                bdata.NJPB_DT = Tmp{6}; clear uBT Tmp;
            end
        end
    end 
    pp.bdata = bdata;
end