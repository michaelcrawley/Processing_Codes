function [varargout] = NF_PostProcessing_AtP(trigger,varargin)

    %This function will reorder signals in actuation-blocks into
    %processing-blocks. The program assumes that all the signals to be
    %reshaped are organized in the same manner. Inputs need to be in the
    %top level, i.e. nf (rather than nf.ablocks.smp). This program will
    %automatically work on .smp only (since .p makes no sense).
    
    %Get Organization Constants
    S = size(varargin{1}.pblocks.p);
    NPp = S(1);
    NV = length(varargin);
    NB = length(trigger.a_i); %number of raw processing blocks
    [NPa,~,~] = size(varargin{1}.ablocks.smp); %size of actuation points in subblocks and number of channels
    RtS = S(3)/NB; %raw to sub-block conversion
    nNBps = floor(NPa*length(trigger.a_i{1})/NPp)*NB; %number of smoothed processing sub-blocks
    StS = nNBps/NB; %Smoothed to sub-block conversion
    
    %Initialize outputs
    varargout = varargin;
    
    %Replace blocks
    for q = 1:NV
        [~,NC,~] = size(varargin{q}.ablocks.smp);
        varargout{q}.pblocks.smp = zeros(NPp,NC,nNBps);
        ablocknum = 1;
        pblocknum = 1;
        for n = 1:NB
            sig = varargout{q}.pblocks.p(:,:,RtS*(n-1)+(1:RtS));
            sig = permute(sig,[1 3 2]);
            sig = reshape(sig,NPp*RtS,[]);
            for k = 1:length(trigger.a_i{n})
                sig(trigger.a_i{n}(k)+(0:NPa-1),:) = varargin{q}.ablocks.smp(:,:,ablocknum);
                ablocknum = ablocknum +1;
            end
            gb = rem(NPp*RtS-trigger.a_i{n}(1),NPp)+1;
            sig = sig(trigger.a_i{n}(1):end-gb,:);
            sig = reshape(sig,NPp,[],NC);
            varargout{q}.pblocks.smp(:,:,pblocknum+(0:StS-1)) = permute(sig,[1 3 2]);
            pblocknum = pblocknum+StS;
        end
    end
end