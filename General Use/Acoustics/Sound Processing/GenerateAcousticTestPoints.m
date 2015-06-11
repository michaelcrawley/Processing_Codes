function [] = GenerateAcousticTestPoints(Stdf, m, PW, IA, AO,R)
% Generates text file containing test points for use with automatic control
% VI for acoustic DAQ.
% Inputs:
%   Stdf:   1D array of Strouhal numbers for forcing
%   m:      Nx2 matrix of forcing modes (use NaN in second column for
%           single mode forcing)
%   PW:     Pulse width in us (use -1 or 'auto' for auto dutycycle)
%   IA:     Increment actuators (optional)
%   AO:     Actuators On (user needs to convert from bin to dec) (optional)
%   R:      Rotate actuators (optional)

    N = 8; %number of channels for control system
    filename = [date,' testpoints.txt'];
    [A,B] = size(m);
    if (A ~= 2) && (B ~= 2)
        error('Mode matrix incorrectly defined');
    elseif (A == 2) && (B ~= 2)
        m = m';
        [A,~] = size(m);
    end    
    if strcmpi(PW,'auto') 
        PW = -1;
    end
    if ~exist('IA','var')
        IA = 0;
    end
    if ~exist('AO','var')
        AO = (2^N)-1;
    elseif strcmpi(AO,'all')
        AO = (2^N)-1;
    end
    if ~exist('R','var')
        R = 0;
    end
    
    testpoints = zeros(length(Stdf)*A*length(PW)*length(IA)*length(AO)*length(R),7);
    place = 1;
    for i = 1:length(R)
       for j = 1:length(AO)
          for k = 1:length(IA)
            for l = 1:length(PW)
               for n = 1:A
                   for o = 1:length(Stdf)
                      testpoints(place,:) = [Stdf(o) m(n,:) PW(l) IA(k) AO(j) R(i)];
                      place = place+1;
                   end
               end
            end              
          end
       end        
    end 
    if exist(filename,'file') ~= 2 %add baselines at beginning and end of test
        testpoints = [ zeros(1,7); testpoints; zeros(1,7) ]; 
    else
        testpoints = [ testpoints; zeros(1,7) ];
    end
    fid = fopen(filename,'a');
    fprintf(fid,'%6.2f \t %i \t %i \t %i \t %i \t %i \t %i \n',testpoints');
    fclose(fid);
end