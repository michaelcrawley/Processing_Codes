function data = grabSpectrum_v2(multi,varargin)
% data = grabSpectrum_v2(multi);
% data = grabSpectrum_v2(multi,SD);
%
%Choose single grab or multiple grab
%   multi = 'off' == single grab
%   multi = 'on'  == multiple grab

if isempty(varargin)
   SD = cd;
else
    SD = varargin{1};
end
[FileName,PathName] = uigetfile('*.fft*','Select Spectrum File',SD,'MultiSelect',multi);

if iscell(FileName)
    for n = 1:length(FileName)
        data{n} = dlmread([PathName FileName{n}],'\t', 1, 0);
        disp(['Grabbing: ' FileName{n}])
    end
else
    data = dlmread([PathName FileName],'\t', 1, 0);
    disp(['Grabbing: ' FileName])
end