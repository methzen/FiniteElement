function solveur(TOL,itermax,FORCEXTG,SDISPT,PROP,XYZ,LE)

global DISPTD
[NUMNP, NDOF] = size(XYZ);
 NEQ = NDOF*NUMNP+size(LE,1);
%
iter=0;
RESIDU=FORCEXTG;
conv=max(abs(RESIDU))/(max(abs(FORCEXTG))+1);
%			
fprintf('\n iter conv');
while conv>TOL && iter<itermax
     [GKF,~]=HYPER3D0(PROP, 1, XYZ,DISPTD,LE);
     SOL=solve(GKF,RESIDU,SDISPT,NEQ,iter,PROP);
     DISPTD=DISPTD+SOL;
     [~,FORCE]=HYPER3D0(PROP, 0, XYZ,DISPTD,LE);
     RESIDU=FORCEXTG-FORCE;
     conv=checkconvergence(RESIDU,SDISPT,NEQ,FORCEXTG);
     iter=iter+1 ;
     fprintf('\n %3d  %12.3e',iter,conv);
end
%% affichage
NDOF=3*size(XYZ,1);
DSP=full(DISPTD(1:NDOF));
DSP=reshape(DSP,3,size(XYZ,1));
XYZ2=XYZ+10*DSP';
for i=1:size(LE,1)
plotHexa8(XYZ(LE(i,:),:),'b')
plotHexa8(XYZ2(LE(i,:),:),'r')
end

end