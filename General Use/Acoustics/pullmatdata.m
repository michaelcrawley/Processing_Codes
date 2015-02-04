[vars vardir] = uigetfile('Multiselect','on');
for j = 1:length(vars)
   load([vardir '\' vars{j}]);
   vars{j} = vars{j}(1:end-4);
end
N = length(vars);
comp.dOASPL = zeros(N,length(eval([vars{1},'.processing_params.ChPol'])));
comp.dOASPL_DT = comp.dOASPL;
comp.dAAE = zeros(1,N);
comp.dAAE_DT = comp.dAAE;
comp.dNJP = comp.dAAE;
comp.dNJP_DT = comp.dAAE;
comp.ff = comp.dAAE;
comp.Stdf = comp.dAAE;
comp.ChPol = eval([vars{1},'.processing_params.ChPol']);

comp.f = eval([vars{1},'.data.F_axis']);
comp.Std = eval([vars{1},'.data.Std_axis']);
for j = 1:N
    comp.ff(j) = eval([vars{j},'.data.ff']);
    comp.Stdf(j) = eval([vars{j},'.data.Stdf']);
    comp.dOASPL(j,:) = eval([vars{j},'.results.dOASPL']);
    comp.dOASPL_DT(j,:) = eval([vars{j},'.results.dOASPL_DT']);
    comp.dAAE(j) = eval([vars{j},'.results.dAAE']);
    comp.dAAE_DT(j) = eval([vars{j},'.results.dAAE_DT']);
    comp.dNJP(j) = eval([vars{j},'.results.dNJP']);
    comp.dNJP_DT(j) = eval([vars{j},'.results.dNJP_DT']);    
end


figure;
contourf(comp.ChPol,comp.Stdf,comp.dOASPL_DT);colorbar;xlabel('Polar Angle');ylabel('St_d_f');colormap cool;

figure;
plot(comp.Stdf,comp.dAAE_DT);xlabel('St_d_f');ylabel('\DeltaAAE (dB)');

figure;
plot(comp.Stdf,comp.dNJP_DT);xlabel('St_d_f');ylabel('\DeltaNJP (dB)');

clear N j vardir vars;