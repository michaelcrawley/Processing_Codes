function [varargout] = NF_PostProcessing_PtA(trigger,varargin)
    %Note that this function is coded to only work with .smp blocks!
    %This function will reorder signals in process-blocks into
    %actuation-blocks. The program assumes that all the signals to be
    %reshaped are organized in the same manner. Inputs need to be in the
    %bottom level, i.e. nf.pblocks.smp
    
    %Get organizational constants
    NV = length(varargin);
    NB = length(trigger.a_i); %number of data sampling blocks
    NPa = round(mean(cell2mat(cellfun(@(x) diff(x),trigger.a_i,'uniformoutput',false)'))); %number of points in actuation period
    [NPp,~,~] = size(varargin{1}); %size of processing points in subblocks and number of channels
    nNBp = floor((NPa-1)*length(trigger.a_i{1})/NPp); %number of processing sub-blocks
    nNBa = floor(NPp*nNBp/NPa); %number of actuation sub-blocks
    NBt = nNBa*NB; %total number of actuation blocks
    
    %Reshape blocks    
    for q = 1:NV
        [~,NC,~] = size(varargin{q});
        varargout{q} = zeros(NPa,NC,NBt);
        ablocknum = 1;
        for n = 1:NB
            pblocknum = (n-1)*nNBp+(1:nNBp);
            a_i = trigger.a_i{n}-trigger.a_i{n}(1)+1;
            a_i = a_i(a_i <= nNBp*NPp-NPa);
            sig = varargin{q}(:,:,pblocknum);
            sig = permute(sig,[1 3 2]);
            sig = reshape(sig,[],NC);
            for k = 1:length(a_i)
                varargout{q}(:,:,ablocknum) = sig(a_i(k)+(0:NPa-1),:);
                ablocknum = ablocknum +1;
            end
        end
        varargout{q} = varargout{q}(:,:,1:ablocknum-1);
    end
end