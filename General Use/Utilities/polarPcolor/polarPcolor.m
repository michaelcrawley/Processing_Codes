function polarPcolor(R,theta,Z) 

%%
%                           ------------------
%                           ------------------
%                                  GOAL
%                           ------------------
%                           ------------------
% Similar to pcolor but in polar coordinates

%                           ------------------
%                           ------------------
%                                  INPUT
%                           ------------------
%                           ------------------ 

% R : type: float
%     size: [1 x Nrr ] where Nrr = numel(R)
%     dimension: radial distance (m)

% theta : type: float
%       size: [1 x Ntheta ] where Ntheta = numel(theta)
%       dimension: azimuth or elevation angle (deg)
%       The zero is defined with respect to the North

% Z : type: float
%     size: [Ntheta x Nrr]
%     dimension: user's defined

%                           ------------------
%                           ------------------
%                                  OUTPUT
%                           ------------------
%                           ------------------

% none: a figure is displayed

%%
% Exemple 1

% R = linspace(25,50,100);
% theta = linspace(5,75,50)
% Z = rand(50,100);
% polarPcolor(R,theta,Z)

% Exemple 2
% R = linspace(3,10,100);
% theta = linspace(0,180,360);
% Z = linspace(0,10,360)'*linspace(0,10,100);
% 
% polarPcolor(R,theta,Z)
%%
%                           ------------------
%                           ------------------
%                                 SCRIPT
%                           ------------------
%                           ------------------


%%                          ------------------
%                           ------------------
%                             Dimension check
%                           ------------------
%                           ------------------

Nrr = numel(R);
Ntheta = numel(theta);
test = size(Z);
if ~isequal(test,[Ntheta,Nrr]),
    Z=Z';
    test = size(Z);
    if ~isequal(test,[Ntheta,Nrr]),
        error(' dimension of Z does not agree with dimension of R and Theta')
    end
end


%%
%                           ------------------
%                           ------------------
%                             MESH DEFINITION
%                           ------------------
%                           ------------------
Nangle = 10; % number of displayed graduation for the angle
Nradius = 7; % number of displayed graduation for the radius

% Definition of the mesh 
radiusMesh = linspace(min(R),max(R)*1.01,Nradius);
angleMesh = linspace(min(theta),max(theta),Nangle);
Radius_range = radiusMesh(end) - radiusMesh(1); % get the range for the radius  
Radius_normalized = R/Radius_range; %[0,1]

% get hold state
cax = newplot;  
ls = get(cax,'gridlinestyle');
hold on;

% define the grid in polar coordinates
angleGrid = linspace(90-min(theta),90-max(theta),100);
xGrid = cosd(angleGrid);
yGrid = sind(angleGrid); 

% transform data in polar coordinates to Cartesian coordinates.
YY = (Radius_normalized)'*cosd(theta);
XX = (Radius_normalized)'*sind(theta);

% plot data on top of grid 
pcolor(XX,YY,Z','parent',cax);
shading flat
      
set(gca,'dataaspectratio',[1 1 1]);axis off;  

                       
%%
%                           ------------------
%                           ------------------
%                              DRAW CIRCLES
%                           ------------------
%                           ------------------




    %------------------------------------
    % Plot the graduation for the radius
    %------------------------------------
    contour = abs((radiusMesh - radiusMesh(1))/Radius_range+R(1)/Radius_range); 
%   whos contour xGrid yGrid
    
    % plot circles
    for kk=1:length(contour) 
        plot(xGrid*contour(kk), yGrid*contour(kk),ls,'color','black','linewidth',1);
    end
     
      %plot radius
        cost = cosd(90-angleMesh); % the zero angle is aligned with North
        sint = sind(90-angleMesh); % the zero angle is aligned with North   
        for kk = 1:length(angleMesh)
            plot(cost(kk)*contour,sint(kk)*contour,ls,'color','black','linewidth',1,...
             'handlevisibility','off');
        % plot graduations of angles
        text(1.07.*contour(end).*cost(kk),1.07.*contour(end).*sint(kk),...
            sprintf('%.3g^{o}',angleMesh(kk)),...
            'horiz', 'center', 'vert', 'middle');
        end    
    
    angleRange= abs(max(theta)-min(theta));
        for kk=1:Nradius 
            text((contour(kk)).*cosd(90-mean(angleMesh)),(contour(kk)).*sind(90-mean(angleMesh)),...
            num2str(radiusMesh(kk),2),'verticalalignment','bottom',...
        'handlevisibility','off','parent',cax);
        end         
        
        text(contour(end).*1.1.*cosd(90-mean(angleMesh)),contour(end).*1.1.*sind(90-mean(angleMesh)),...
        ['distance (m)'],'verticalalignment','bottom',...
        'handlevisibility','off','parent',cax);
%     end        
     
    % colorbar  
    caxis([quantile(Z(:),0.05),quantile(Z(:),0.95)])
    c =colorbar;
    axes(c) 
	location = 'NorthOutside';
	set(c,'fontSize',18,'location',location,'Fontweight','bold');
	title(' your variable here ')  

      



    % set figure to screen size
    set(gcf,'Position',get(0,'ScreenSize'))
    set(gcf,'color','w');

    % Increase the size of the text on the figure
    angleGrid = findall(gcf,'Type','text');
    for ii = 1:length(angleGrid),
        set(angleGrid(ii),'FontSize',18,'Fontweight','bold')

    end  
    
end

