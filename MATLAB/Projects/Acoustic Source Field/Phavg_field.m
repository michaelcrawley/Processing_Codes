function [varargout] = Phavg_field(lafpa,varargin)

    indx = [find(diff(lafpa)>2.25)+1; length(lafpa)];
    L = max(diff(indx));
    
    varargout = cell(length(varargin),1);
    for q = 1:length(varargin)
        [M,N,~] = size(varargin{q});
        varargout{q} = zeros(M,N,L);
        num_average = zeros(L,1);
        
        for k = 1:(length(indx)-1)
            varargout{q}(:,:,1:(indx(k+1)-indx(k))) = varargout{q}(:,:,1:(indx(k+1)-indx(k))) + varargin{q}(:,:,indx(k):(indx(k+1)-1));
            num_average(1:(indx(k+1)-indx(k))) =  num_average(1:(indx(k+1)-indx(k))) +1;
        end
        
        for k = 1:L
            varargout{q}(:,:,k) = varargout{q}(:,:,k)/num_average(k);
        end
    end

end