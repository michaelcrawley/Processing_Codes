function [tx,sx,c] = waveletAnalysis(d,LF,FS,wvtype)


% LF = 100;   %Hz - lowest desired frequency
% FS = 200000;    %Hz - signal sampling frequency
% wvtype = 'mexh'; %Wavelet type
% d = dlmread([cd '/Work/m09_baseline_t30.3.nos'],'\t'); d = d(1:8192,9);

s = (2:FS/LF);
c = cwt(d,s,wvtype);

tx = (0:1/FS:(length(d)-1)/FS);
sx = FS./s;

