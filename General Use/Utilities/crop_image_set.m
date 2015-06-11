function crop_image_set(src_dir,flist, varargin)

    if nargin == 2
        prefix = '';
        suffix = '_1';
        [~,~,ext] = fileparts(flist{1});
        cmds = '';
    elseif nargin == 3
        prefix = varargin{1};
        suffix = '';
        [~,~,ext] = fileparts(flist{1}); 
        cmds = '';
    elseif nargin == 4
        prefix = varargin{1};
        suffix = varargin{2};
        [~,~,ext] = fileparts(flist{1});
    elseif nargin == 5
        prefix = varargin{1};
        suffix = varargin{2};
        ext = varargin{3};
    end
    
    %read in first image, have user identify cropping area
    img = imread([src_dir filesep flist{1}]);
    [subdir,filename] = fileparts(flist{1});
    h = figure; 
    imshow(img);
%     jframe = get(h,'JavaFrame');
%     jframe.setMaximized(true);
    [x,y] = ginput(2);
    close(h);
    
    %create index arrays
    if x(1) < x(2)
        X = floor(x(1)):ceil(x(2));
    else
        X = floor(x(2)):ceil(x(1));
    end
    
    if y(1) < y(2)
        Y = floor(y(1)):ceil(y(2));
    else
        Y = floor(y(2)):ceil(y(1));
    end
   
    %save first image
    cimg = img(Y,X,:);
    imwrite(cimg,[src_dir filesep subdir filesep prefix filename suffix ext]);
    
    for n = 2:length(flist)
        img = imread([src_dir filesep flist{n}]);
        [subdir,filename] = fileparts(flist{n});
        cimg = img(Y,X,:);
        imwrite(cimg,[src_dir filesep subdir filesep prefix filename suffix ext]);
    end
end