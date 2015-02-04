function [ average, lengthscale, N,ensembleaverage, ReynoldsStress, autocor, crosscor, Spectra] = H2( data )
    %Completed by Michael Crawley 04/26/2010 for ME 813


    dt = mean(diff(data.t));
    %%Part 1
    figure
    subplot(4,1,1),plot(data.t,data.pwall/1000);xlabel('t (s)');ylabel('pwall (kPa)');
    subplot(4,1,2),plot(data.t,data.p1/1000);xlabel('t (s)');ylabel('p1 (kPa)');
    subplot(4,1,3),plot(data.t,data.u1);xlabel('t (s)');ylabel('u1 (m/s)');
    subplot(4,1,4),plot(data.t,data.v1);xlabel('t (s)');ylabel('v1 (m/s)');
    
    %%Part 2
    
    average.pwall = mean(data.pwall);
    average.p1 = mean(data.p1);
    average.u1 = mean(data.u1);
    average.v1 = mean(data.v1);
    
    %%Part 3
    [autocor.pwall s.a.pwall] = xcorr(data.pwall-average.pwall,'coeff');
    [autocor.p1 s.a.p1] = xcorr(data.p1-average.p1,'coeff');
    [autocor.u1 s.a.u1] = xcorr(data.u1-average.u1,'coeff');
    [autocor.v1 s.a.v1] = xcorr(data.v1-average.v1,'coeff');
    l = length(autocor.pwall(1:(end/2)))+1;
    
    I.pwall = find(autocor.pwall(l:end) < 0,1)+l-2;
    I.p1 = find(autocor.p1(l:end) < 0,1)+l-2;
    I.u1 = find(autocor.u1(l:end) < 0,1)+l-2;
    I.v1 = find(autocor.v1(l:end) < 0,1)+l-2;
    
    lengthscale.pwall = dt*trapz(autocor.pwall(l:I.pwall));
    lengthscale.p1 = dt*trapz(autocor.p1(l:I.p1));
    lengthscale.u1 = dt*trapz(autocor.u1(l:I.u1));
    lengthscale.v1 = dt*trapz(autocor.v1(l:I.v1));
    
    figure
    subplot(4,1,1),plot(dt*s.a.pwall(l:I.pwall),autocor.pwall(l:I.pwall), [ 0 lengthscale.pwall lengthscale.pwall], [1 1 0]);xlabel('t (s)');title('Pwall');ylabel('Autocorrelation Coefficient');ylim([ 0 1.1]);
    subplot(4,1,2),plot(dt*s.a.p1(l:I.p1),autocor.p1(l:I.p1), [ 0 lengthscale.p1 lengthscale.p1], [1 1 0]);xlabel('t (s)');title('P1');ylabel('Autocorrelation Coefficient');ylim([ 0 1.1]);
    subplot(4,1,3),plot(dt*s.a.u1(l:I.u1),autocor.u1(l:I.u1), [ 0 lengthscale.u1 lengthscale.u1], [1 1 0]);xlabel('t (s)');title('U1');ylabel('Autocorrelation Coefficient');ylim([0 1.1]);
    subplot(4,1,4),plot(dt*s.a.v1(l:I.v1),autocor.v1(l:I.v1), [ 0 lengthscale.v1 lengthscale.v1], [1 1 0]);xlabel('t (s)');title('V1');ylabel('Autocorrelation Coefficient');ylim([0 1.1]);
    
    %%Part 4
    N = floor(data.t(end)/lengthscale.u1);
    
    %%Part 5
    n = length(data.t);
    dx = round(linspace(1,n,N));
    ensembleaverage.u1 = mean(data.u1(dx));
    ensembleaverage.v1 = mean(data.v1(dx));
    ensembleaverage.pwall = mean(data.pwall(dx));
    ensembleaverage.p1 = mean(data.p1(dx));
    ReynoldsStress.uu = mean((data.u1(dx)-ensembleaverage.u1).^2);
    ReynoldsStress.vv = mean((data.v1(dx)-ensembleaverage.v1).^2);
    ReynoldsStress.uv = mean((data.u1(dx)-ensembleaverage.u1).*(data.v1(dx)-ensembleaverage.v1));
    
    %%Part 6
    figure
    plot(data.u1(dx)-ensembleaverage.u1,data.v1(dx)-ensembleaverage.v1,'*');xlabel('u''');ylabel('v''');
    
    %%Part 7
    [crosscor.p1u1 s.c.p1u1] = xcorr(data.p1-average.p1,data.u1-average.u1,'coeff');
    [crosscor.p1v1 s.c.p1v1] = xcorr(data.p1-average.p1,data.v1-average.v1,'coeff');
    [crosscor.p1pwall s.c.p1pwall] = xcorr(data.p1-average.p1,data.pwall-average.pwall,'coeff');
    figure
    subplot(3,1,1),plot(dt*s.c.p1u1,crosscor.p1u1);xlabel('t (s)'); ylabel('Correlation Coefficient');title('Crosscorrelation between P1 and U1');xlim([-.005 .005]);
    subplot(3,1,2),plot(dt*s.c.p1v1,crosscor.p1v1);xlabel('t (s)'); ylabel('Correlation Coefficient');title('Crosscorrelation between P1 and V1');xlim([-.005 .005]);
    subplot(3,1,3),plot(dt*s.c.p1pwall,crosscor.p1pwall);xlabel('t (s)'); ylabel('Correlation Coefficient');title('Crosscorrelation between P1 and Pwall');xlim([-.005 .005]);
    
    %%Part 8
    edt = mean(diff(data.t(dx)));
    NFFT = 2^nextpow2(length(dx));
    omega = 0.5*(1/edt)*linspace(0,1,NFFT/2+1);
    Spectra.pwall = 2*abs(fft(xcorr(data.pwall(dx) - ensembleaverage.pwall),NFFT));
    Spectra.p1 = 2*abs(fft(xcorr(data.p1(dx) - ensembleaverage.p1),NFFT));
    Spectra.u1 = 2*abs(fft(xcorr(data.u1(dx) - ensembleaverage.u1),NFFT));
    Spectra.v1 = 2*abs(fft(xcorr(data.v1(dx) - ensembleaverage.v1),NFFT));
    figure
    subplot(4,1,1),loglog(omega,Spectra.pwall(1:NFFT/2+1));title('Pwall Spectrum'); xlabel('\omega');ylabel('E');
    subplot(4,1,2),loglog(omega,Spectra.p1(1:NFFT/2+1));title('P1 Spectrum');xlabel('\omega');ylabel('E');
    subplot(4,1,3),loglog(omega,Spectra.u1(1:NFFT/2+1));title('U1 Spectrum');xlabel('\omega');ylabel('E');
    subplot(4,1,4),loglog(omega,Spectra.v1(1:NFFT/2+1));title('V1 Spectrum');xlabel('\omega');ylabel('E');
    
end