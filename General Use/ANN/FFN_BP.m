function [nneFun, MSE, weights, phi, params] = FFN_BP(xt,d,arch,varargin)
%This code serves as a engine for a static feed-forward,
%back-propagation artificial neural network. The network is defined to be
%fully connected (this cannot be changed by user).
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
%               'mode':     Processing mode, either 'batch' or 'pattern'.
%                           Default is 'batch'.
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
    dw = cell(L-1,1); %weight updates
    dw_old = dw; %storage for old weight updates (used for momentum in learning)
    y = cell(L,1); %Container for activation levels, training inputs will be in first level
    v = cell(L-1,1); %Container for neuronal summation levels
    delta = cell(L-1,1);   
    
    %Get Parameters
    [maxepoch,eta,alpha,econv_total,econv_change,mode,phi,dphi,weights] = get_options(arch,varargin);
    params = struct('maxepoch',maxepoch,'eta',eta,'alpha',alpha,'econv_total',econv_total,'econv_change',econv_change','mode',mode,'phi',phi,'dphi',dphi,'weights',weights);
    
    %Set processing parameters
    if strcmpi(mode,'batch')
        Npass = 1;
        Nset = N;
        flag = true;
    elseif strcmpi(mode,'pattern')
        Npass = N;
        Nset = 1;
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
        %Randomize training set order
        if ~flag
            rorder = randperm(N);
            xt = xt(rorder,:);
            d = d(rorder,:);
        end

        for q = 1:Npass
            %Feed-Forward Computations
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            if flag
                y{1} = xt;
            else
                y{1} = xt(q,:);
            end
            input = [y{1} -ones(Nset,1)]; %concat bias term
            for nn = 1:L-1
                %compute neuron summation
                v{nn} = input*weights{nn};
                %compute activation level
                y{nn+1} = phi(v{nn});
                %concat input to account for neuron bias
                input = [y{nn+1} -ones(Nset,1)];
            end

            %Back-Propagation Computations
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %Compute Mean-squared error
            if flag
                e = d-y{end}; %error signal at output level
            else
                e = d(q,:) - y{end};
            end
            MSE(n) = MSE(n) + mean(sum(e.^2,2))/2/Npass; 

            %Compute Output-Layer delta
            delta{end} = e.*dphi(v{end});
            if n == 1
                dw{end} = eta*(delta{end}.'*[y{end-1} -ones(Nset,1)]).'/Nset; %include bias weight updating
            else
                dw{end} = eta*(delta{end}.'*[y{end-1} -ones(Nset,1)]).'/Nset+alpha*dw_old{end};
            end

            %Compute Hidden-Layer delta
            for nn = (L-1):-1:2
                e = delta{nn}*weights{nn}(1:end-1,:).';
                delta{nn-1} = e.*dphi(v{nn-1});
                if n == 1
                    dw{nn-1} = eta*[y{nn-1} -ones(Nset,1)].'*delta{nn-1}/Nset;
                else
                    dw{nn-1} = eta*[y{nn-1} -ones(Nset,1)].'*delta{nn-1}/Nset+alpha*dw_old{nn-1};
                end
            end

            %Update weights
            for nn = 1:L-1
                weights{nn} = weights{nn}+dw{nn};
                dw_old{nn} = dw{nn};
            end 
        end

        %%Convergence Computations
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        if (MSE(n) < econv_total) && (abs((MSE(n-1)-MSE(n))/MSE(n-1)) < econv_change)
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
    str = 'phi([x,-ones(size(x,1),1)]*weights{1})';
    for n = 2:L-1
        str = ['phi([' str ',-ones(size(x,1),1)]*weights{',num2str(n),'})'];
    end
    nneFun = eval(['@(x) ' str]);
end

function [maxepoch,eta,alpha,econv_total,econv_change,mode,phi,dphi,weights] = get_options(arch,commands)
    %Determine inputs
    options = {'slope','amplitude','maxepoch','lrp','momentum','wrange','etotal','echange','mode','phi','weights'};
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
    
    %Processing mode
    if any(cmd == 9)
        loc = find(cmd == 9)+1;
        mode = commands{loc};
    else
        mode = 'batch';
    end

    %Activation function
    if any(cmd == 10)
        loc = find(cmd == 10);
        phi = commands{loc+1};
        dphi = commands{loc+2};
    else
        phi = @(v) A*tanh(s*v); %scales between -A and A
        dphi = @(v) A*s*sech(s*v).^2; %derivative with respect to v
    end
    
    %Initial weights
    if any(cmd == 11)
        loc = find(cmd == 11)+1;
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