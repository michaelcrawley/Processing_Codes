function data = grabSpectrum_v2(varargin)

if isempty(varargin)
   SD = cd;
else
    SD = varargin{1};
end
[FileName,PathName] = uigetfile('*.fftNOS','Select Spectrum File',SD,'MultiSelect','on');
for n = 1:length(FileName)
    data{n} = dlmread([PathName FileName{n}],'\t', 1, 0);
    disp(['Grabbing: ' FileName{n}])
end
