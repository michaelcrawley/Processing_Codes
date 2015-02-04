function [flist_out Uout x y] = PreProcessPIV2d(src_dir,flist_in)

    L = length(flist_in); 
    q = 3;

    %read in first file to get sizing information, initialize variables
    data = readimx([src_dir filesep flist_in{1}]);
    [x,y,u,v] = showimx(data);
    
    %throw out garbage data
    x = x(q:end,1);
    y = y(1,:);
    u = u(q:end,:);
    v = v(q:end,:);
    
    [M N] = size(u);
    U = zeros(2*M*N,L); 
    U(:,1) = [u(:); v(:)];
    
    %read in data from all files
    for n = 2:L
        data = readimx([src_dir filesep flist_in{n}]);
        [~,~,u,v] = showimx(data); 
        u = u(q:end,:);
        v = v(q:end,:);
        U(:,n) = [u(:); v(:)];       
    end
    
    %Identify any duplicate images
    R = corrcoef(U)-eye(L); %find correlations between images, subtract out autocorrelation values
    I = []; %index array for garbage images
    for n = 1:L
        if all(n ~= I)
            new = find(R(:,n) >= 0.99); %find images identical to image 'n'
            I = [I; new]; %add garbage image indices to array
        end
    end
    
    %Identify any blank (or nearly blank) images
    Usum = sum(abs(U));
    new = find(Usum < 100);
    I = [I; new];
    
    %Identify any random vector data sets based on correlation to mean
    imgs = true(1,L);
    for n = 1:length(I)
        if any(I(n) == (1:L));
            imgs(I(n)) = 0;
        end
    end
    umean = mean(U(:,imgs),2);
    X = [U umean];
    R = corrcoef(X);
    correlation = R(:,end);
    new = find(correlation < 0.3); %find images that do not correlate well with mean velocity
    I = [I; new];
    
    for n = 1:length(new)
        if any(new(n) == (1:L));
            imgs(new(n)) = 0;
        end
    end
    Uout = U(:,imgs);

    %Rename garbage images (changes extension from VC7 to garbage)
    for n = 1:length(I)
        [~,name] = fileparts([src_dir filesep flist_in{I(n)}]);
        movefile([src_dir filesep flist_in{I(n)}],[src_dir filesep name '.garbage']);
    end
    flist_out = flist_in(imgs);
end