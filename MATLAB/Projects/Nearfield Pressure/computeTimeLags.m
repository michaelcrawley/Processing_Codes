function [tau_ac tau_con d] = computeTimeLags(X,Y,xc,a,conv,R,phi)
%This function computes the acoustic and hydrodynamic/acoustic propagation
%time lags from nearfield points to the specified far-field microphone.
%Make sure units are consistent!

    %Compute distance to far-field microphone
    r = [R*cosd(phi) R*sind(phi)]; %X,Y position of far-field mic
    d = sqrt((r(1)-X).^2+(r(2)-Y).^2); %distance from near-field to far-field mics
    
    %Compute acoustic propagation path time lag
    tau_ac = d/a; %acoustic raypath
    
    %Compute time lag from jet centerline at end of potential core to
    %observer
    tau_const = sqrt((r(1)-xc).^2+r(2).^2)/a;
    
    %Compute hydrodynamic/acoustic propagation path time lag
    [Mc,Nc] = size(conv);
    [M,N] = size(X);
    c = repmat(conv,M/Mc,N/Nc);
    tau_con = tau_const+(xc-X)./c; %convective-acoustic raypath
end