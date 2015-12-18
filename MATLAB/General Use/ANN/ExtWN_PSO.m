function [nneFun, MSE, w] = ExtWN_PSO(xt,d,arch,varargin)
%This code constructs a static, fully-connected wavelet network. 
%Learning is accomplished solely using particle swarm
%optimization. The user must open a worker pool prior to calling this
%function if parallelization is desired.
%Last updated on 2015-04-21 by Michael Crawley.
%Inputs:
%           xt:         [M,N] matrix of training sets, where M is the
%                       number of training examples and N is the number of 
%                       input neurons.
%           d:          [M,P] matrix of desired outputs for training, where
%                       M is the number of training sets and P is the
%                       number of output neurons.
%           arch:       Number of hidden wavelons.
%           options:
%               'maxepoch':     Maximum number of epochs for training. Default
%                               is 1e3.
%               'wrange':       Range for random assignent for initial weights.
%                               Default is +/- 2.4/N * 10.
%               'etotal':       MSE error for convergence. Default is eps^1/2.
%               'phi':          Activation function. Default is hyperbolic tangent.
%               'weights':      Initial weights, default is random within 
%                               [-wrange,wrange].
%               'particles':    Number of particles in swarm. Default is 1e3
%               'velocities':   Velocity update parameters - See
%                               constructor function.
%Outputs:
%           nneFun:     Function that performs neural network computations
%                       provided a single input set.
%           MSE:        Mean squared error over epochs.
%           weights:    Final neural network weights.
%           phi:        Activation function used by neural network.
%           params:     Struct containing all network parameters

    %Get matrix information
    [N,a1] = size(xt); %get number of training sets
    [N2,a2] = size(d);
    if N ~= N2
        error('Input Matrix Mismatch');
    end
    arch = [a1;arch(:);a2]; %append input and output neurons
    
    %Initialize Parameters
    [wavelet,phi,maxepoch,econv_total,swarm,np,chi,alpha,c] = get_options(arch,varargin);
    MSE = zeros(maxepoch,1); %Container for Mean-Squared-Error over epochs  
    L = length(arch);
    
    %Build neural net function   
    nneFun = build_neural_function(arch,wavelet,phi,N);
    
    %Initialize local and global best values and positions
    for n = 1:np
        swarm(n).value = mean(mean(nneFun(xt,swarm(n).position,d).^2,2),1);
        swarm(n).best_value = swarm(n).value;
        swarm(n).best_position = swarm(n).position;
    end
    [global_best.value,I] = min([swarm.best_value]);
    global_best.position = swarm(I).best_position;
    
    %Set up progress bar
    cpb = ConsoleProgressBar();
    cpb.setLeftMargin(4);   % progress bar left margin
    cpb.setTopMargin(1);    % rows margin
    cpb.setLength(40);      % progress bar length: [.....]
    cpb.setMinimum(0);      % minimum value of progress range [min max]
    cpb.setMaximum(100);    % maximum value of progress range [min max]
    cpb.setPercentPosition('left');
    cpb.setTextPosition('right');
    cpb.setMinimum(0);
    cpb.setMaximum(maxepoch);
    cpb.start();
    
    for n = 1:maxepoch %epoch counter 
        
        %%FeedForward Computations
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        for q = 1:np
            %Update Position and Velocity
            for k = 1:L-1
                swarm(q).velocity{k} = chi*(alpha*swarm(q).velocity{k} + c(1)*rand(1)*(swarm(q).best_position{k}-swarm(q).position{k}) + c(2)*rand(1)*(global_best.position{k}-swarm(q).position{k}));
                swarm(q).position{k} = swarm(q).position{k} + swarm(q).velocity{k};
            end
            
            %Update error and local best value & position
            swarm(q).value = mean(mean(nneFun(xt,swarm(q).position,d).^2,2),1);
            if swarm(q).value < swarm(q).best_value
                swarm(q).best_value = swarm(q).value;
                swarm(q).best_position = swarm(q).position;
            end
        end
                
        %%Update MSE and Global bests        
        [itr_best_value,I] = min([swarm.best_value]);
        if itr_best_value < global_best.value
            global_best.value = itr_best_value;
            global_best.position = swarm(I).position;
        end
        MSE(n) = itr_best_value;
        
        %%Convergence Computations
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        if (MSE(n) < econv_total)
            MSE = MSE(1:n); %truncate
                                
            text = sprintf('Iteration: %d/%d [Converged]', n, maxepoch);  
            cpb.setValue(n);  	% update progress value
            cpb.setText(text);  % update user text
            fprintf('\n');
            break;
        end
        
        %Update user if not converged
        if n > 1
            dmse = MSE(n)-MSE(n-1);
        else
            dmse = 0;
        end
        text = sprintf('Iteration: %d/%d [MSE:%.3g,change:%.3g]', n, maxepoch,MSE(n),dmse);
        cpb.setValue(n);  	% update progress value
        cpb.setText(text);  % update user text
    end
    fprintf('\n');  

    %Build Final Function
    w = global_best.position; 
    nneFun = build_final_neural_function(arch,wavelet,phi,w);
end

function [nne_out] = build_neural_function(arch,wavelet,phi,N)
    L = length(arch);
    wavelons = @(x,w) prod(wavelet((permute(repmat(x,[1,1,arch(2)]),[1 3 2])-repmat(permute(w{1}(:,:,2),[3 2 1]),[N,1,1]))./repmat(permute(w{1}(:,:,1),[3 2 1]),[N,1,1])),3);
    nne_out = @(x,w) phi([wavelons(x,w),x,-ones(N,1)]*w{2});
    for n = 3:L-1
        nne_out = @(x,w) phi([nne_out(x,w),-ones(N,1)]*w{n});
    end
    nne_out = @(x,w,d) nne_out(x,w) - d;
end

function [nne_out] = build_final_neural_function(arch,wavelet,phi,w)
    L = length(arch);
    wavelons = @(x) prod(wavelet((permute(repmat(x,[1,1,arch(2)]),[1 3 2])-repmat(permute(w{1}(:,:,2),[3 2 1]),[size(x,1),1,1]))./repmat(permute(w{1}(:,:,1),[3 2 1]),[size(x,1),1,1])),3);
    nne_out = @(x) phi([wavelons(x),x,-ones(size(x,1),1)]*w{2});
    for n = 3:L-1
        nne_out = @(x) phi([nne_out(x),-ones(size(x,1))]*w{n});
    end 
end

function [wavelet,phi,maxepoch,econv_total,swarm,np,chi,alpha,c] = get_options(arch,commands)
    %Determine inputs
    options = {'slope','amplitude','maxepoch','wrange','etotal','phi','particles','weights','velocity'};
    cmd = zeros(1,length(commands));
    for n = 1:length(commands)
        if ischar(commands{n})
            cmd(n) = find(ismember(options,commands{n}));
        end
    end
    
    %slope for activation function
    if any(cmd == 1)
        loc = find(cmd == 1)+1;
        if isa(commands{loc},'function_handle')
            wavelet = commands{loc};
            dwavelet = commands{loc+1};
        else
            mother = commands{loc};
            m = commands{loc+1};
            [wavelet,dwavelet] = MotherWavelets(mother,m);
        end 
    else
        mother = 'mexican';
        [wavelet,dwavelet] = MotherWavelets(mother);
    end
    
    %amplitude for activation function
    if any(cmd == 2)
        loc = find(cmd == 2)+1;
        A = commands{loc};
    else
        A = 1.0;
    end
    
    %maximum number of epochs for training
    if any(cmd == 3)
        loc = find(cmd == 3)+1;
        maxepoch = commands{loc};
    else
        maxepoch = 1e3; 
    end
    
    %range for random weight initialization
    if any(cmd == 4)
        loc = find(cmd == 4)+1;
        wrange = commands{loc};
    else
        wrange = 4.8/max(arch(1),1)*10; %fix in case no hidden layers exist
    end
    
    %Error limit for convergence based off of MSE
    if any(cmd == 5)
        loc = find(cmd == 5)+1;
        econv_total = commands{loc};
    else
        econv_total = sqrt(eps);
    end
    
    %Activation function
    if any(cmd == 6)
        loc = find(cmd == 6);
        phi = commands{loc+1};
    else
        phi = sidestep_memory_phi(A); %scales between -A and A
    end
    
    %Number of particles in swarm
    if any(cmd == 7)
        loc = find(cmd == 7)+1;
        np = commands{loc};
    else
        np = 1e3;
    end
    
    %Initial weights
    if any(cmd == 8)
        loc = find(cmd == 8)+1;
        swarm = commands{loc};
    else
        L = length(arch);
        swarm(np) = struct('position',[],'velocity',[],'value',[],'best_position',[],'best_value',[]); %initialize
        parfor q = 1:np
            swarm(q).position = cell(L-1,1); %use different cell for each layer
            swarm(q).position{1}(:,:,1) = rand(arch(1),arch(2)) - 0.5;
            swarm(q).position{1}(:,:,2) = rand(arch(1),arch(2)) - 0.5;
            swarm(q).position{2} = rand(arch(2)+arch(1)+1,arch(3)) - 0.5;         
            swarm(q).velocity = cell(L-1,1); %use different cell for each layer
            swarm(q).velocity{1}(:,:,1) = rand(arch(1),arch(2)) - 0.5;
            swarm(q).velocity{1}(:,:,2) = rand(arch(1),arch(2)) - 0.5;
            swarm(q).velocity{2} = rand(arch(2)+arch(1)+1,arch(3)) - 0.5;
            
            for n = 3:L-1
                %weights are randomly initialized to be between (-wrange,wrange)
                %velocities are initialized as being double this
                %column length is incremented by one for bias term 
                swarm(q).position{n} = wrange*rand(arch(n)+1,arch(n+1))-wrange/2; 
                swarm(q).velocity{n} = wrange*rand(arch(n)+1,arch(n+1))-wrange/2; 
            end
        end
    end
    
    %Swarm velocity parameters
    if any(cmd == 9)
        loc = find(cmd == 8)+1;
        chi = commands{loc}(1);
        alpha  = commands{loc}(2);
        c = commands{loc}(3:4);
    else
        chi = 0.98;
        alpha = 0.75;
        c = [0.125 0.125];
    end

end

function out = sidestep_memory_phi(A)
    out = @(v) A*tanh(v);
end

function [wavelet,dwavelet] = MotherWavelets(mother,m)
%This function returns a function handle, 'wavelet', given a dimension and
%string for the name of the mother wavelet. A parameter to modify the base
%mother wavelet can also be provided. For standard mother wavelets, the
%dimension provided will be an integer, for Spatio-Temporal mother
%wavelets, the dimension will be a string beginning with 'ST' and
%including the total number of dimensions.

    if ~exist('m','var'), m = []; end
    switch lower(mother)
        case 'morlet'
            if isempty(m), m = 6; end
            wavelet = @(x) exp(1i*m*x).*exp(-0.5*x^2);
            dwavelet = @(x) 1i*m*exp(1i*m*x).*exp(-0.5*x^2) - x.*exp(1i*m*x).*exp(-0.5*x^2);
        case 'mexican'
            wavelet = @(x) (1-x.^2).*exp(-0.5*x.^2);
            dwavelet = @(x) (x.^3 - 3*x).*exp(-0.5*x.^2);
        otherwise
            error('Undefined Mother Wavelet');
    end
end
    