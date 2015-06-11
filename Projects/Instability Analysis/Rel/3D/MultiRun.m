function MultiRun(m,nmode,omg,alph,Rt,tr,workers,sleepstate)
    %Used to run multiple parameter sets of AxiJet using parallel
    %processing toolbox
    %Version 1.3
    !powercfg -change -standby-timeout-ac 0
    savedir = uigetdir(pwd,'Specify save folder');
    
    if matlabpool('size') > 0
        matlabpool close
    end
    matlabpool(workers);
    
    parfor i = 1:length(m)
        disp(['Processing for Mach(',num2str(m(i),2),') AzMode(',num2str(nmode(i)),') Rt(',num2str(Rt(i)),') Tr(',num2str(tr(i),2),')']);
        AxiJetTv1(m(i),nmode(i),omg{1,i},alph(i),Rt(i),tr(i),savedir);
        disp(['AxiJet completed for Mach(',num2str(m(i),2),') AzMode(',num2str(nmode(i)),') Rt(',num2str(Rt(i)),') Tr(',num2str(tr(i),2),')']);
    end
    sleep = ['powercfg -change -standby-timeout-ac ',num2str(sleepstate)];
    dos(sleep);
    matlabpool close
end