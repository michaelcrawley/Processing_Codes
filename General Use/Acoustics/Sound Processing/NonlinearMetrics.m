function [results] = NonlinearMetrics(rRaw,results,pp)
    
    %reshape pressure waveform data
    rRaw = rRaw';
    rRaw = rRaw(:);
    rRaw = reshape(rRaw,pp.Nch,pp.BS,[]); %each channel is now a row vector, and the block number is the third dimension

    h =  1/pp.sampleRate;
    coefm = mNumericalDerivative(1,2,h,pp.BS);
    [~,~,N] = size(rRaw);
    
    %initialize variables
    Pskewness = zeros(pp.Nch,N);
    Pkurtosis = Pskewness; dPskewness = Pskewness; dPkurtosis = Pskewness;
    
    for n = 1:N
        drRaw = rRaw(:,:,n)*coefm; %calculate first temporal derivative of pressure signal (in voltage)

        %calculate metrics in voltage domain
        Pskewness(:,n) = skewness(rRaw(:,:,n),1,2);
        Pkurtosis(:,n) = kurtosis(rRaw(:,:,n),1,2);
        dPskewness(:,n) = skewness(drRaw,1,2);
        dPkurtosis(:,n) = kurtosis(drRaw,1,2);
    end

    %average metrics for all blocks
    results.Pskewness = mean(Pskewness,2);
    results.Pkurtosis = mean(Pkurtosis,2);
    results.dPskewness = mean(dPskewness,2);
    results.dPkurtosis = mean(dPkurtosis,2);
end