function InterpSpatialArray(x,y,t,trig,sig,dx,y1,optset,cpus)

    %Set default parameters - these will likely not change between the
    %cases
    phys.Uj = 285.99; 
    phys.D = 0.0254;
    phys.M = 0.9;
    phys.Temp = 251.31; %static jet temperature
    phys.theta = 8.6; %angle of the array with respect to jet axis
    phys.To = phys.Temp*(1 + .2*phys.M^2);
    phys.a = sqrt(1.4*287*phys.To);
    phys.FS = 1/(mean(diff(t)));
    %far-field locations to grab - 30^o & 90^o (roughly)
    phys.ff_xi = [length(x),1];
    phys.ff_yi = [length(y),length(y)];
    phys.xLES = x;
    phys.yLES = y;
    matversion = 1.21;
    if exist('optset','var') && ~isempty(optset)
        fields = fieldnames(optset);
        for n = 1:length(fields)
            phys.(fields{n}) = optset.(fields{n});
        end
    end
    
    %format trigger
    trigger.sig = trig;
    if ~isempty(trig)
        ai = find(diff(trig) > 0) +1;
        if (ai(1) - unique(diff(ai))) > 0
            ai = [1; ai];
        end
        trigger.a_i{1} = ai;
    end
       
    %create desired grid
    phys.x = min(x(:)):dx:max(x(:)); 
    phys.y = tand(phys.theta)*(phys.x-1) + y1;

    %Open Matlab Pool if Necessary
    if ~exist('cpus','var')||isempty(cpus), cpus = 0; end
    poolobj = gcp('nocreate');
    if ~isempty(poolobj) && cpus ~= poolobj.NumWorkers
        delete(poolobj);
        poolobj = parpool(cpus);
    elseif isempty(poolobj) && cpus > 0
        poolobj = parpool(cpus);
    end

    N = length(t);
    intwvfm = zeros(N,length(phys.x));
    parfor n = 1:N
        intwvfm(n,:) = interp2(x,y,sig(:,:,n),phys.x,phys.y,'spline');
    end    
    nf.lblocks.p = intwvfm;
    
    ff.lblocks.p = zeros(N,length(phys.ff_xi));
    for n = 1:length(phys.ff_xi)
        ff.lblocks.p(:,n) = squeeze(sig(phys.ff_yi(n),phys.ff_xi(n),:));
    end
  
    fname = ['M',num2str(phys.M),'_r',num2str(y1),'_a',num2str(phys.theta),'.mat'];
    save(fname,'nf','ff','phys','trigger','matversion');
    delete(poolobj);
end