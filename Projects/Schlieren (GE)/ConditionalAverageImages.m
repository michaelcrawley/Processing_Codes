function[R] = ConditionalAverageImages(src_dir,flist,cHi,cWi,dthreshold)

    if ~exist('alpha','var'), dthreshold = 0.9; end    
    files = cellfun(@(file) fullfile(src_dir,file),flist,'UniformOutput',0);
   
    M = length(cHi);
    N = length(cWi);
    nF = length(flist);
    
    chk = zeros(M*N,length(flist));
    for n = 1:nF
        img = double(imread(files{n}));
        chk(:,n) = reshape(img(cHi,cWi),M*N,1);
    end
    
    R = corrcoef(chk); %correlation coefficient
    
    done = false; 
    LL = 0.2;
    threshold = dthreshold;
    nImg = nF/10; %define required number of images for correlation as at least 10% of set
    while ~done
        alpha = R > threshold;
        L = sum(alpha);   %Finds the number of images with correlation above "threshold"
        if max(L) > nImg   %If the necessary number of images is found
            done = true;
            CorrCoeff = threshold;
            disp([num2str(max(L)) ' Positive Match Found at Corr. Coeff. = ' num2str(threshold)])
        elseif threshold < LL   %If threshold has fallen below acceptable lower limit
            done = true;
            CorrCoeff = 0;
            disp('No Positive Result')
        else
            threshold = threshold-0.01;
        end
    end
    
    [~,I] = max(L);
    matches0 = alpha(:,I);
    Ph0 = flist(matches0);
    
%     done = false;
%     threshold = -dthreshold;
%     while ~done
%         alpha = R(:,I) < threshold;
%         L = sum(alpha);
%         if L > N
%             done = true;
%             CorrCoeff(2) = threshold;
%             disp([num2str(L) ' Negative Match Found at Corr. Coeff. = ' num2str(threshold)])
%         elseif threshold > -LL
%             done = true;
%             CorrCoeff(2) = 0;
%             disp('No Positive Result')
%         else
%             threshold = threshold + 0.01;
%         end
%     end
    
%     matches180 = alpha;
%     Ph180 = flist(matches180);

    AverageImages(src_dir,Ph0,['CAvg Ph0,nImg',num2str(length(Ph0))]);
%     AverageImages(src_dir,Ph180,['CAvg Ph180,nImg',num2str(length(Ph180))]);
end