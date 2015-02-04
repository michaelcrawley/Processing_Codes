%This program will take the exported vector field of a baseline, upstream
%case and read the data, calculate the freestream velocity and all of the
%boundary layer properties and print them to the screen.
clear
%close all
clc

global r
global To
global Po
global u_tau


%%%%%% THE VALUES BELOW THIS POINT NEED TO BE DEFINED %%%%%%


%set sw to 0 to automatically determine the Mach number based on the given
%stagnation temperature, or set sw to 1 to use the predetermined Mach
%number
sw=1;

%constants to define
M_infty=2.3;%1.89; %Freestream Mach number
xshift=0; %(+) shifts the zero to the right
yshift=0; %(+) shifts the zero up
To=298.15; %Stagnation temperature (K)
Po=48;%41; %psig
u_tau=16.79; %Friction velocity m/s lookup
lm=20; %The spanwise +/- limit over which to average for the u-profile (mm)
rect=[-10,10,20,10]; %The area from which to calculate the freestream velocity/Mach#,
                    %[left,right,top,bottom] (mm)
percent=98; %Percent of freestream to limit boundary layer thickness
r=0.896;%0.949; %compressibility factor for the calculation of momentum and
                %displacement thickness
rm=1; %First point from the wall to keep
xloc=94.49; %Streamwise location from the end of the nozzle (mm)
h=2.87;%1.5; %tunnel height at the nozzle exit (in)


%%%%%%% THE VALUES ABOVE THIS POINT NEED TO BE DEFINED %%%%%%

Po=(Po+14.696)/14.696*101325; %converting to Pa

% Velocity data from file
i=1;
fprintf('Requesting user input.\n');
[d fs] = userSelectFiles('Select the Vector Fields to be Averaged','.mat');

% Get current file
fi = fullfile(d, fs{i});
load(fi);
global dir_name
dir_name=d;

global logf
logf=fopen(fullfile(d,'BL_data.txt'),'w');

[yp,zp]=meshgrid(y,z);
xp=x;
clear x y z
z{i}=zp;
y{i}=yp;
x{i}=xp;
u{i}=Up';
v{i}=Vp';
W{i}=Wp';
c{i}=sqrt(1.4*287*To*(1+0.2*M_infty^2)^-1);

%[z{i},y{i},x{i},u{i},v{i},w{i},c{i}]=read_vel_data(xshift,yshift);

if sum(size(x{1}))~=4   %Checking to ensure the file was properly read
    %Rectangular area

    [v_infty,c_infty]=mean_vel_rect(z{1},y{1},u{1},c{1},rect(1),rect(2),rect(3),rect(4));
    if ~sw
        M_infty=v_infty/c_infty;
    end
    fprintf('The freestream Mach number is %.2f\n',M_infty);
    fprintf('The freestream velocity is %.1f m/s\n\n',v_infty);
    
    fprintf(logf,'The freestream Mach number is %.2f\n',M_infty);
    fprintf(logf,'The freestream velocity is %.1f m/s\n\n',v_infty);

    %Horizontally RMS Averaged Velocity
    [v_max{i}, v_prof{i}]=horizontal_vrms(z{i},y{i},u{i},lm);

    %Calculating the boundary layer properties
    wd=size(y{1});
    Mcorr=h*25.4/tan(asin(1/M_infty));
    xloc=xloc+Mcorr;
    delta=bl_prop(v_infty,M_infty,xloc,percent,[0,v_prof{1}(rm+1:size(v_prof{1},2))],[0,y{1}(1,rm:wd(2)-1)]);
    mmBLplot(delta,y{i},M_infty,v_prof{i},v_infty,u_tau);
    % vprof3d(u{1},y{1},z{1},delta,v_infty);
    % run temp_prof_plotter
    fclose(logf);
end