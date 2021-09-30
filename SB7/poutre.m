%
% Two-element example page 131
%
clear all
clc
format long g
%% Donnee d'entrees
% Nodal coordinates
%geom=Plaque(2,0.5,0.1,3,1,1);
%XYZ=full(geom.Nodes3D(:,2:4));
% % %
% % % Element connectivity
%LE=geom.Connect3D;
LE=[1 2 3 4 5 6 7 8];
%
 XYZ=[ 0 0 0;
       2 0 0;
       2 2 0;
       0 2 0;
       0 0 2;
       2 0 2;
       2 2 2;
       0 2 2];
%  LE=[1 2 3 4 5 6 7 8 ] ;
% %   B=gradtransf(XYZ,[0 0 0])
for i=1:size(LE,1)
plotHexa8(XYZ(LE(i,:),:),'b')
end
%AfficheRepere(0,0,0,eye(3),0)
%% Material properties HYPERELASTIQUE
%
E=1;
v=0.0;
% LAM=E*v/((1+v)*(1-2*v));
% MU=E/(2+2*v);
% a10=MU/2;
PROP=[E   v];
%% Conditions limites
% External forces [Node, DOF, Value]
F=0.25/4;
EXTFORCE=[2 1 F
          3 1 F
          6 1 F
          7 1 F];

% Prescribed displacements [Node, DOF, Value]

SDISPT=[1	1	0
1	2	0
1	3	0
4	1	0
4   2   0
4	3	0
5	1	0
5	2	0
5   3   0
8	1	0
8   2   0
8   3   0
];
%% Assemblage de la force exterieure
[NUMNP, NDOF] = size(XYZ);				
NE = size(LE,1);
NEQ = NDOF*NUMNP;
FORCEXTG = sparse(NEQ+NE,1); % FORCE EXTERIEUR GENERALISEE
DISPTD= sparse(NEQ+NE,1);
if size(EXTFORCE,1)>0
LOC = NDOF*(EXTFORCE(:,1)-1)+EXTFORCE(:,2);
FORCEXTG(LOC) = FORCEXTG(LOC) + EXTFORCE(:,3);
end
%
TOL=1.0e-6;
itermax=5; 
iter=0;
RESIDU=FORCEXTG;
conv=max(RESIDU)/(max(FORCEXTG)+1.0);
%
%[GKF,~]=HYPER3D(PROP, 1, XYZ,DISPTD,LE);

while conv>TOL && iter<itermax
     [GKF,~]=HYPER3D0(PROP, 1, XYZ,DISPTD,LE);
     SOL=solve(GKF,RESIDU,SDISPT,NEQ+NE,iter);
     DISPTD=DISPTD+SOL;
     [~,FORCE]=HYPER3D0(PROP, 1, XYZ,DISPTD,LE);
     RESIDU=FORCEXTG-FORCE;
     conv=checkconvergence(RESIDU,SDISPT,NEQ+NE,FORCEXTG)
     iter=iter+1     
end

NDOF=3*size(XYZ,1);
DSP=full(DISPTD(1:NDOF));
DSP=reshape(DSP,3,size(XYZ,1));
XYZ2=XYZ+10*DSP';
for i=1:size(LE,1)
plotHexa8(XYZ(LE(i,:),:),'b')
%[~, SHPD, ~] = SHAPEL([0 0 0], XYZ2(LE(i,:),:));
%F=DSP(:,LE(i,:))*SHPD' + eye(3)
%E=(F'*F-eye(3))/2
plotHexa8(XYZ2(LE(i,:),:),'r')
end