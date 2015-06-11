%   This code is to verify results in Prof. Koochesfahani's paper
%   incom-huerre.m --- inviscid incompressible disturbance solution
%   search for complex eigenvalue alpha for specified range of real omg
%==========summary of computation process==================================
%   There are a range of real omg(s)
%   For every omg, search for the right alpha using prediction-and-try method
%   For every guessed alpha, use RKF45 meathod to integrate from y=y1 to y=y2
%==========summary of computation process==================================

clear
close all
pause(1)

% data ================================================
% frequency range (omg)
%omgall=0.2:-0.01:0;mode 2
omgall=0.1%:0.001:0.78; %mode 1

% omgall = 0.52;
% alph = complex(0.2,-.3);

num=length(omgall);
% initial guessed alph
%alph=complex(0.7,-0.09); mode 2
alph=complex(0.2,-0.09); %mode 1
% computation domain
y1=-8;
y2=8;
% error controls
nmax=200;
tol=1.e-5;
abserr=1.e-5;
% data for specific omg to be save for plot
nomg=1;
PYN=0;          % PYN=1 if data found for this omg
% data ================================================

% Loop of omg ############################################
% Matrics to save true alph and coresponding omg
talph=zeros(1,num);
tomg=zeros(1,num);
nn=0;  % number of omg whose right alph is found

multiWaitbar('Frequency Range:',0,'Color',[0.1 0.5 0.8]);
multiWaitbar('Iteration:',0,'Color',[0.2 0.9 0.3]);

for kk=1:num
    omg=omgall(kk);
    disp(['omg=',num2str(omg),', process=',num2str(kk),'/',num2str(num)])
    % extrapolate guessed starting alph for this omg
    if kk>2
        galphr=interp1(tomg(1:nn),real(talph(1:nn)),omg,'spline','extrap');
        galphi=interp1(tomg(1:nn),imag(talph(1:nn)),omg,'spline','extrap');
        alph=complex(galphr,galphi);
    end
    % Look for alph for this omg
    multiWaitbar('Iteration:','Reset');
    multiWaitbar('Iteration:','Color',[0.2 0.9 0.3]);
    for itr=1:nmax+1 
        multiWaitbar('Iteration:','Increment',1/(nmax+1));
        % use guessed alph for c
        c=omg/alph;
        % Impose Dirichlet BC at y=y1
        z1=zeros(4,1);
        z1(1)= exp(alph*y1);
        z1(2)= alph*z1(1);
        z1(3)= y1*z1(1);
        z1(4)= alph*z1(3)+z1(1);
        % intergrate from y1 to y2, output (z,y) at every integration point
        [zz,yy]=RKF45(z1,y1,y2,abserr,alph,c);
        % check accuracy of z at y=y2, guess new alph
        er=alph*zz(1,end)+zz(2,end);
        der=zz(1,end)+alph*zz(3,end)+zz(4,end);
        dalph=-er/der;
        alpherr=abs(dalph)/abs(alph);      
        disp(['itr=',' ',num2str(itr),', ','error=',num2str(alpherr),', abs error=',num2str(abs(zz(1,end)))])
        if itr>10
           multiWaitbar('Iteration:','Color',[0.8 0.0 0.1]); 
        end
        if (alpherr<tol)
            % right alph found, end search
            multiWaitbar('Iteration:','Value',1);
            pause(1/60);
            break            
        else
            % new alph for next trial
            alph=alph+dalph;
        end        
    end
    multiWaitbar('Frequency Range:','Increment',1/num);


    % right alph found or not
    if itr>nmax
        % right alph not found, break loop of omg
        disp(['search failed: No alph found for this omg in ',num2str(nmax),' trials']);disp(' ')
        break
    else
        % right alph found
        nn=kk;
        disp(['alph for current omg is: ', num2str(alph)])
        disp(['# of trials to find alph for this omg= ', num2str(itr)]); disp(' ')
        talph(kk)=alph;
        tomg(kk)=omg;
        % For specific omg, save z (if found)
        if kk==nomg
            PYN=1;
            omgsave=omg;
            alphsave=alph;
            zsave=zz(1,:);
            ysave=yy;
        end
    end
end
multiWaitbar('CLOSEALL');
% Loop of omg ############################################



% plot alph against omg
if nn>0
    % alph found for at least one omg
    % plot alphr against omg
    alphr=real(talph);
    alphi=imag(talph);
    figure(1)
    plot(tomg(1:nn),alphr(1:nn));
    axis([0,1,0,2])
    title('\alphar for a range of \omega')
    xlabel('\omega')
    ylabel('\alphar')
    grid on
    % plot alphi against omg
    figure(2)
    plot(tomg(1:nn),-alphi(1:nn));
    axis([0,1,0,0.4])
    title('-\alphai for a range of \omega')
    xlabel('\omega')
    ylabel('-\alphai')
    grid on    
else
    disp('No alph found for any omg')
end



% plot eigenfunction for the specific omg
if PYN==1
    % plot eigenfunction (complex form)
    figure(3)
    plot(zsave)
    title(['Eigenfunction (disturbance v) when \omg= ',num2str(omgsave)])
    xlabel('\phir') ; ylabel('\phii')
    grid on    
    % plot real part of eigenfunction against y
    figure(4)
    zreal=real(zsave);
    plot(zreal,ysave)
    title('real part of eigenfunction')
    xlabel('\phir'); ylabel('y');
    grid on    
    % plot imaginary part of eigenfunction against y
    figure(5)
    zimag=imag(zsave);
    plot(zimag,ysave)
    title('imaginary part of eigenfunction')
    xlabel('\phii'); ylabel('y');
    grid on    
    disp(' ')
    disp(['alph for omg [',num2str(omgsave),'] is [',num2str(alphsave),']']);disp(' ')
else
    % alph for this specific omg not found
    disp(' ')
    omgspec=omgall(nomg);
    disp(['No alph found for the specific omg [',num2str(omgspec),'], No data to plot']);disp(' ')
end




