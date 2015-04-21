function [nneFun, MSE] = WN_BP(xt,d,arch,varargin)
%This code implements a multidimensional wavelet network per Alexandridis
%2013 (Neural Networks). Learning is performed using the standard
%backpropagation algorithm in batch mode. The network is fully connected, 
%and includes linear connections from the input to the output. This code
%restricts the user to a single populated hidden layer.
%Last updated on 2015-04-14 by Michael Crawley.
%[UNFINISHED]
%Inputs:
%           xt:         [M,N] matrix of training sets, where M is the
%                       number of training examples and N is the number of 
%                       input neurons.
%           d:          [M,P] matrix of desired outputs for training, where
%                       M is the number of training sets and P is the
%                       number of output neurons.
%           arch:       [1,L] array describing network architecture, where
%                       L is the number of layers (excluding input and
%                       output). For a network with three inputs, one output,
%                       and two hidden layers with two neurons, the network
%                       architecture would be [2,2].
%           options:
%               'mother':   Mother wavelet, default is Morlet-6
%               'A':        Amplitude for the default activation function.
%                           Default is 1.0.
%               'maxepoch': Maximum number of epochs for training. Default
%                           is 1e6.
%               'lrp':      Learning-rate-parameter. Default is 0.1.
%               'momentum': Momentum parameter. Default is 0.0.
%               'etotal':   MSE error for convergence. Default is eps^1/2.
%               'weights':  Initial weights, default is random within 
%                           [-wrange,wrange].
%Outputs:
%           nneFun:     Function that performs neural network computations
%                       provided a single input set.
%           MSE:        Mean squared error over epochs.
%           weights:    Final neural network weights.
%           phi:        Activation function used by neural network.
%           params:     Struct containing all network parameters

    %Get training matrix information
    [N,a1] = size(xt); %get number of training sets
    [N2,a2] = size(d);
    if N ~= N2
        error('Input Matrix Mismatch');
    end
    
    %Initialize variables
    arch = [a1;arch(:);a2]; %append input and output neurons
    [wavelet,dwavelet,maxepoch,eta,alpha,econv_total,weights] = get_options(arch,varargin);
    weights_old = weights; %initialize here so that momentum can be included immediately
    
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

    %Perform Back-Propagation Algorithm
    MSE = zeros(maxepoch,1); %Container for Mean-Squared-Error over epochs
    for n = 1:maxepoch %epoch counter 
        
        %Feed-Forward Computations
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        
        dilation = repmat(permute(weights{1},[3 2 1]),[N,1,1]);
        translation = repmat(permute(weights{2},[3 2 1]),[N,1,1]);
        inputs = repmat(xt,[1,1,arch(2)]);
        z_ij = (inputs-translation)./dilation;
        
        wavelets = wavelet(z_ij);
        wavelons = prod(wavelets,3); %single level of hidden wavelons only
        output = wavelons*weights{3} + xt*weights{4} + weights{5}; %includes linear terms from input and bias terms in addition to the outputs of the wavelons
        

        %Back-Propagation Computations
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %Compute Mean-squared error
        e = d-output; %error signal at output level
        MSE(n) = mean(sum(e.^2,2))/2; 

        %Compute partial derivatives of error terms
        de{5} = mean(e,1); %bias terms
        de{4} = xt'*e/N; %pass-through terms
        de{3} = wavelons'*e/N; %wavelon-output terms
        de{2} = -weights{3}./weights{1} uhhhh; %wavelon scale terms
        de{1} = (z_ij.*de{2}).'*e/N; %wavelon translation terms            

        %Update weights
        swap = weights;
        for n = 1:5
            weights{n} = weights{n} + eta*de{n} + alpha*(weights{n}-weights_old{n});
        end
        weights_old = swap;
            

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
    nneFun = @(x) prod(wavelet((repmat(x,arch(1),1)-weights{2})./weights{1}),2).'*weights{3} + x*weights{4} + weights{5};
end


function [wavelet,dwavelet,maxepoch,eta,alpha,econv_total,weights] = get_options(arch,commands)
    %Determine inputs
    options = {'mother','maxepoch','lrp','momentum','etotal','weights'};
    cmd = zeros(1,length(commands));
    for n = 1:length(commands)
        if ischar(commands{n})
            cmd(n) = find(ismember(options,commands{n}));
        end
    end
    
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
    
    
    %maximum number of epochs for training
    if any(cmd == 2)
        loc = find(cmd == 3)+1;
        maxepoch = commands{loc};
    else
        maxepoch = 1e6; 
    end

    %learning-rate-parameter
    if any(cmd == 3)
        loc = find(cmd == 4)+1;
        eta = commands{loc};
    else
        eta = 0.1;
    end
    
    %momentum parameter
    if any(cmd == 4)
        loc = find(cmd == 5)+1;
        alpha = commands{loc};
    else
        alpha = 0.0;
    end
    
    %Error limit for convergence based off of MSE
    if any(cmd == 5)
        loc = find(cmd == 6)+1;
        econv_total = commands{loc};
    else
        econv_total = sqrt(eps);
    end
    
    %Initial weights
    if any(cmd == 6)
        loc = find(cmd == 7)+1;
        weights = commands{loc};
    else
        %Weights will organized in five levels:
        %   Level 1: scale weight for hidden wavelons
        %   Level 2: translation weight for hidden wavelons
        %   Level 3: Linear amplitude weight for hidden-output neurons
        %   Level 4: Linear amplitude weight for input-output
        %   Level 5: Bias weights at output
        weights{1} = rand(arch(1),arch(2)) - 0.5;
        weights{2} = rand(arch(1),arch(2)) - 0.5;
        weights{3} = rand(arch(2),arch(end)) - 0.5;
        weights{4} = rand(arch(1),arch(end)) - 0.5;
        weights{5} = rand(arch(end),1)-0.5;
    end
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