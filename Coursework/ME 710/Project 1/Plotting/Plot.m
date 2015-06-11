%%Function for creating plots

%T1
X = a;
h = figure;
contourf(100*x,100*x',Tn(:,:,X)); colorbar; colormap gray;xlabel('x (cm)');ylabel('y (cm)');title(['Tnumeric (K) for time ',num2str(round(t(X))),' (s)']);
saveas(h,'Tn t1','fig'); saveas(h,'Tn t1','png');close(h);
h = figure;
contourf(100*x,100*x',Ta(:,:,X)); colorbar; colormap gray;xlabel('x (cm)');ylabel('y (cm)');title(['Tanalytic (K) for time ',num2str(round(t(X))),' (s)']);
saveas(h,'Ta t1','fig');saveas(h,'Ta t1','png'); close(h);
h = figure;
contourf(100*x,100*x',Ta(:,:,X)-Tn(:,:,X)); colorbar; colormap gray;xlabel('x (cm)');ylabel('y (cm)');title(['Error (K) for time ',num2str(round(t(X))),' (s)']);
saveas(h,'Terror t1','fig'); saveas(h,'Terror t1','png');close(h);

%T2
X = b;
h = figure;
contourf(100*x,100*x',Tn(:,:,X)); colorbar; colormap gray;xlabel('x (cm)');ylabel('y (cm)');title(['Tnumeric (K) for time ',num2str(round(t(X))),' (s)']);
saveas(h,'Tn t2','fig'); saveas(h,'Tn t2','png');close(h);
h = figure;
contourf(100*x,100*x',Ta(:,:,X)); colorbar; colormap gray;xlabel('x (cm)');ylabel('y (cm)');title(['Tanalytic (K) for time ',num2str(round(t(X))),' (s)']);
saveas(h,'Ta t2','fig'); saveas(h,'Ta t2','png');close(h);
h = figure;
contourf(100*x,100*x',Ta(:,:,X)-Tn(:,:,X)); colorbar; colormap gray;xlabel('x (cm)');ylabel('y (cm)');title(['Error (K) for time ',num2str(round(t(X))),' (s)']);
saveas(h,'Terror t2','fig'); saveas(h,'Terror t2','png');close(h);

%T3
X = c;
h = figure;
contourf(100*x,100*x',Tn(:,:,X)); colorbar; colormap gray;xlabel('x (cm)');ylabel('y (cm)');title(['Tnumeric (K) for time ',num2str(round(t(X))),' (s)']);
saveas(h,'Tn t3','fig');saveas(h,'Tn t3','png'); close(h);
h = figure;
contourf(100*x,100*x',Ta(:,:,X)); colorbar; colormap gray;xlabel('x (cm)');ylabel('y (cm)');title(['Tanalytic (K) for time ',num2str(round(t(X))),' (s)']);
saveas(h,'Ta t3','fig'); saveas(h,'Ta t3','png'); close(h);
h = figure;
contourf(100*x,100*x',Ta(:,:,X)-Tn(:,:,X)); colorbar; colormap gray;xlabel('x (cm)');ylabel('y (cm)');title(['Error (K) for time ',num2str(round(t(X))),' (s)']);
saveas(h,'Terror t3','fig'); saveas(h,'Terror t3','png');close(h);

%T4
X = d;
h = figure;
contourf(100*x,100*x',Tn(:,:,X)); colorbar; colormap gray;xlabel('x (cm)');ylabel('y (cm)');title(['Tnumeric (K) for time ',num2str(round(t(X))),' (s)']);
saveas(h,'Tn t4','fig'); saveas(h,'Tn t4','png');close(h);
h = figure;
contourf(100*x,100*x',Ta(:,:,X)); colorbar; colormap gray;xlabel('x (cm)');ylabel('y (cm)');title(['Tanalytic (K) for time ',num2str(round(t(X))),' (s)']);
saveas(h,'Ta t4','fig'); saveas(h,'Ta t4','png');close(h);
h = figure;
contourf(100*x,100*x',Ta(:,:,X)-Tn(:,:,X)); colorbar; colormap gray;xlabel('x (cm)');ylabel('y (cm)');title(['Error (K) for time ',num2str(round(t(X))),' (s)']);
saveas(h,'Terror t4','fig'); saveas(h,'Terror t4','png');close(h);

%T5
X = e;
h = figure;
contourf(100*x,100*x',Tn(:,:,X)); colorbar; colormap gray;xlabel('x (cm)');ylabel('y (cm)');title(['Tnumeric (K) for time ',num2str(round(t(X))),' (s)']);
saveas(h,'Tn t5','fig'); saveas(h,'Tn t5','png');close(h);
h = figure;
contourf(100*x,100*x',Ta(:,:,X)); colorbar; colormap gray;xlabel('x (cm)');ylabel('y (cm)');title(['Tanalytic (K) for time ',num2str(round(t(X))),' (s)']);
saveas(h,'Ta t5','fig');saveas(h,'Ta t5','png');  close(h);
h = figure;
contourf(100*x,100*x',Ta(:,:,X)-Tn(:,:,X)); colorbar; colormap gray;xlabel('x (cm)');ylabel('y (cm)');title(['Error (K) for time ',num2str(round(t(X))),' (s)']);
saveas(h,'Terror t5','fig'); saveas(h,'Terror t5','png');close(h);