function saveFigure_v2(FIG,FN,R,varargin)
%
%This program is for automatically saving figures which are of journal
%quality. The chosen image format used in the code is png because of its
%excellent image quality for digitally created images while also employing
%a lossless compression algorithm. TIF produces similar image results, but
%since it has no compression algorithm, the file sizes are easily 10x
%larger. The program crops off unnecessary gray space around figures, sets
%font properties, and saves the image at the desired resolution.
%
%Figures may need to be resized for final submission. This can be
%accomplished one of two ways: the piece of code below (assumes you have 
%the Image Processing Toolbox), or Photoshop. Photoshop produces better
%results and can be done efficiently using macros for batch processing. If
%you wish to use MATLAB resizing, uncomment the code below and put in a
%desired width.
%
%Figures may need to be in TIF format for actual journal submission. If so,
%I'd recommend using Photoshop because it has the capability to save TIF's
%with LZW compression (this capability doesn't exist in MATLAB). As before
%this operation can be easily accomplished using batch processes within
%Photoshop. TIF files with LZW can be opened and handled by many programs
%including publishers.
%
%NOTE: Figure windows should not be shrunk as a method of resizing because
%all objects do not scale by the same amount which can result in
%undesirable overlaps.
%
%GENERAL RULE: Do not manipulate figure with mouse before running this
%program - Undesirable results may occur due to broken relationships
%between figure objects.
%
%
%Program accepts the following input types:
%   saveFigure(FIG,FN,R) - This sets spacing between rows/columns to zero
%   and saves figure to file name FN with resolution R.
%       This is also the proper call for a figure with only one axis.
%   saveFigure(FIG,FN,R,rSpace,cSpace) - Uses specified spacing values.
%
%FIG - handle to desired figure (e.g. gcf)
%FN - file path and name without file extension (png will be appended)
%R (DPI) - figure resolution (600 is recommended for journals)
%cSpace (percentage) - Spacing between columns beyond necessary minimum
%rSpace (percentage) - Spacing between rows beyond necessary minimum
%


    %12 and BOLD produce clearly legible results even after image shrinking
    %based on previous experience. Change if you wish.
FS = 11;        %Font Size
FW = 'bold';    %Font Weight

if ~isempty(varargin)
    rSpace = varargin{1};
    cSpace = varargin{2};
else
    cSpace = 0;
    rSpace = 0;
end

Q = get(FIG,'Children');
if isempty(Q)
    error('No Axis handles found!');    %Spits an error if you try to pass a bad figure handle
end

for n = 1:length(Q) %Resets all Text objects to desired size and weight
    if isempty(strfind(get(Q(n),'Type'),'ui'))
        set(Q(n),'FontSize',FS,'FontWeight',FW);
        set(get(Q(n),'XLabel'),'FontSize',FS,'FontWeight',FW);
        set(get(Q(n),'YLabel'),'FontSize',FS,'FontWeight',FW);
        set(get(Q(n),'ZLabel'),'FontSize',FS,'FontWeight',FW);
        set(get(Q(n),'Title'),'FontSize',FS,'FontWeight',FW);
        set(Q(n),'GridLineStyle','--');

        QQ = findobj(Q(n),'Type','text');
        if ~isempty(QQ)
            for m = 1:length(QQ)
                set(QQ(m),'FontSize',FS-1,'FontWeight',FW); %Uses 1 point smaller font size for text objects to avoid obscuring data
            end
        end
    end
end
drawnow;

% if ~isempty(varargin)   %Removes gray space with appropriate function call
%     removeGraySpace_v2(FIG,rSpace,cSpace); %for subplots
% else
%     removeGraySpace_v2(FIG);   %for single axis
%     
% %     minorGridlineFix(gca);  %FIXES THE GRIDLINE PRINTING ISSUE
% end

% PPM = get(FIG,'PaperPositionMode');
% set(FIG,'PaperPositionMode','auto');    %Sets parameter so that saved image will accurately reflect the on screen dimensions
print('-dpng',['-r' num2str(R)],[FN '.png'])  %recommended resolution is 600 for publications
% set(FIG,'PaperPositionMode',PPM)    %Restores original Paper parameter

%%%% FOR AUTOMATICALLY RESIZING IMAGES %%%%%%%%%%%%%%%%%%%%%%%%%%
% W = 3;  %Image final width (inches)
%
% D = imfinfo([FN '.png']);
% IM = imread([FN '.png']);
% if strcmp(D.ResolutionUnit,'meter')
%     NC = W*0.0254*D.XResolution;
%     IM = imresize(IM,NC/D.Width);
%     imwrite(IM,[FN '.png'],'png');
% else    %Assumes units are in inches
%     NC = W*D.XResolution;
%     IM = imresize(IM,NC/D.Width);
%     imwrite(IM,[FN '.png'],'png');
% end

