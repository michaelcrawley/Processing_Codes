function data = grabSpectrum(varargin)

if isempty(varargin)
   SD = cd;
else
    SD = varargin{1};
end
[FileName,PathName] = uigetfile('*.fftNOS','Select Spectrum File',SD,'MultiSelect','off');
data = dlmread([PathName FileName],'\t', 1, 0);
disp(FileName)