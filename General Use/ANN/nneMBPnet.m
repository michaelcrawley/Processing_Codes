function [nneFun, MSE, weights, phi, params] = nneMBPnet(xt,d,arch,varargin)
%This code serves as a engine for a modified feed-forward,
%back-propagation artificial neural network, which uses the Delta-Bar-Delta 
%learning rule with gradient descent. Per the results shown in Haykin, only
%batch-processing is allowed. 
%Last updated on 2013-03-02 by Michael Crawley.
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
%               's':        Slope for the default activation function. 
%                           Default is 1.0.                 
%               'A':        Amplitude for the default activation function.
%                           Default is 1.0.
%               'maxepoch': Maximum number of epochs for training. Default
%                           is 1e6.
%               'lrp':      Learning-rate-parameter. Default is 0.1.
%               'momentum': Momentum parameter. Default is 0.0.
%               'wrange':   Range for random assignent for initial weights.
%                           Default is +/- 2.4/N.
%               'etotal':   MSE error for convergence. Default is eps^1/2.
%               'echange':  Percent change in MSE for convergence. Default
%                           is 5e-3.
%               'phi':      Activation function; user must give both phi
%                           and dphi/dv. Default is hyperbolic tangent.
%               'weights':  Initial weights, default is random within 
%                           [-wrange,wrange].
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

    %Initialize storage containers
    L = length(arch); 
    Dij = cell(L-1,1); %weight updates
    dw_old = Dij; %storage for old weight updates (used for momentum in learning)
    Sij = Dij;
    Sij_old = Dij;
    y = cell(L,1); %Container for activation levels, training inputs will be in first level
    v = cell(L-1,1); %Container for neuronal summation levels
    delta = cell(L-1,1);   
    
    %Get Parameters
    [maxepoch,eta,alpha,econv_total,econv_change,phi,dphi,weights] = get_options(arch,varargin);
    params = struct('maxepoch',maxepoch,'eta',eta,'alpha',alpha,'econv_total',econv_total,'econv_change',econv_change','phi',phi,'dphi',dphi,'weights',weights);
    
    %Set processing parameters
    if strcmpi(mode,'batch')
        Npass = 1;
        N = N;
        flag = true;
    elseif strcmpi(mode,'pattern')
        Npass = N;
        N = 1;
        flag = false;
    else
        error('Processing mode not set');
    end
    
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

        y{1} = xt;
        input = [y{1} -ones(N,1)]; %concat bias term
        for nn = 1:L-1
            %compute neuron summation
            v{nn} = input*weights{nn};
            %compute activation level
            y{nn+1} = phi(v{nn});
            %concat input to account for neuron bias
            input = [y{nn+1} -ones(N,1)];
        end

        %Back-Propagation Computations
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %Compute Mean-squared error
        e = d-y{end}; %error signal at output level
        MSE(n) = MSE(n) + mean(sum(e.^2,2))/2;         

        %Compute Output-Layer delta
        delta{end} = e.*dphi(v{end});
        Dij{end} = (delta{end}.'*[y{end-1} -ones(N,1)]).'/N; %include bias weight updating

        %Compute Hidden-Layer delta
        for nn = (L-1):-1:2
            e = delta{nn}*weights{nn}(1:end-1,:).';
            delta{nn-1} = e.*dphi(v{nn-1});
            Dij{nn-1} = [y{nn-1} -ones(N,1)].'*delta{nn-1}/N;
        end
        
        %Update learning rate
        if n > 1
            for nn = 1:L-1
                Sij{nn} = (1-xi)*Dij_old{nn} + xi*Sij_old{nn};
                eta{nn} = eta{nn} + kappa.*(Sij_old{nn}.*Dij{nn} > 0) + -beta*eta{nn}.*(Sij_old{nn}.*Dij{nn} < 0);
            end
        end

        %Update weights
        for nn = 1:L-1
            dw_old{nn} = Dij{nn};
            weights{nn} = weights{nn}+eta{nn}.*Dij{nn};            
            Dij_old{nn} = Dij{nn};
        end 


        %%Convergence Computations
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        if (MSE(n) < econv_total) && (abs((MSE(n-1)-MSE(n))/MSE(n-1)) < econv_change)
            MSE = MSE(1:n); %truncate
                                
            text = sprintf('Iteration: %d/%d [Converged]', n, maxepoch);  
            cpb.setValue(n);  	% update progress value
            cpb.setText(text);  % update user text
            sprintf('/n');
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
    sprintf('/n');

    %Build Final Function
    str = 'phi([x,-1]*weights{1})';
    for n = 2:L-1
        str = ['phi([' str ',-1]*weights{',num2str(n),'})'];
    end
    nneFun = eval(['@(x) ' str]);
end

function [maxepoch,eta,alpha,econv_total,econv_change,phi,dphi,weights] = get_options(arch,commands)
    %Determine inputs
    options = {'slope','amplitude','maxepoch','lrp','momentum','wrange','etotal','echange','phi','weights'};
    cmd = zeros(1,length(commands));
    for n = 1:length(commands)
        if ischar(commands{n})
            cmd(n) = find(ismember(options,commands{n}));
        end
    end
    
        %slope for activation function
    if any(cmd == 1)
        loc = find(cmd == 1) + 1;
        s = commands{loc};
    else
        s = 1.0;
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
        maxepoch = 1e6; 
    end

    %learning-rate-parameter
    if any(cmd == 4)
        loc = find(cmd == 4)+1;
        eta = commands{loc};
    else
        eta = 0.1;
    end
    
    %momentum parameter
    if any(cmd == 5)
        loc = find(cmd == 5)+1;
        alpha = commands{loc};
    else
        alpha = 0.0;
    end
    
    %range for random weight initialization
    if any(cmd == 6)
        loc = find(cmd == 6)+1;
        wrange = commands{loc};
    else
        wrange = 4.8/max(arch(1),1); %fix in case no hidden layers exist
    end
    
    %Error limit for convergence based off of MSE
    if any(cmd == 7)
        loc = find(cmd == 7)+1;
        econv_total = commands{loc};
    else
        econv_total = sqrt(eps);
    end
    
    %Error limit for convergence based off of MSE rate of change
    if any(cmd == 8)
        loc = find(cmd == 8)+1;
        econv_change = commands{loc};
    else
        econv_change = 5e-3;
    end
    
    %Activation function
    if any(cmd == 9)
        loc = find(cmd == 9);
        phi = commands{loc+1};
        dphi = commands{loc+2};
    else
        phi = @(v) A*tanh(s*v); %scales between -A and A
        dphi = @(v) A*s*sech(s*v).^2; %derivative with respect to v
    end
    
    %Initial weights
    if any(cmd == 10)
        loc = find(cmd == 10)+1;
        weights = commands{loc};
    else
        L = length(arch);
        weights = cell(L-1,1); %use different cell for each layer
        for n = 1:L-1
            %weights are randomly initialized to be between (-wrange,wrange)
            %column length is incremented by one for bias term 
            weights{n} = wrange*rand(arch(n)+1,arch(n+1))-wrange/2;   
        end
    end
end