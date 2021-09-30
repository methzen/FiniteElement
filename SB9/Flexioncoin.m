%
% Two-element example page 131
%
clear all
clc
format long g
%% Donnee d'entrees
% Nodal coordinates
%geom=Plaque(100,100,1,1,1,1);
%XYZ=full(geom.Nodes3D(:,2:4));
% % %
% % % Element connectivity
%LE=geom.Connect3D;
%LE=[1 2 3 4 5 6 7 8];
%
 XYZ=[ 0 0 0
       100 0 0;
       100 100 0;
       0 100 0;
       0 0 1;
       100 0 1;
       100 100 1;
       0 100 1];
LE=[1 2 3 4 5 6 7 8 ] ;
% %   B=gradtransf(XYZ,[0 0 0])
% for i=1:size(LE,1)
% plotHexa8(XYZ(LE(i,:),:),'b')
% end
%AfficheRepere(0,0,0,eye(3),0)
%% Material properties HYPERELASTIQUE
%
E=200000;
v=0.;
% LAM=E*v/((1+v)*(1-2*v));
% MU=E/(2+2*v);
% a10=MU/2;
PROP=[E   v];
%% Conditions limites
% External forces [Node, DOF, Value]
F=1;
B=LE(1,[1 2 3 4 5 6 7 8]);
EXTFORCE=[3 3 F];

% Prescribed displacements [Node, DOF, Value]
SDISPT=[
B(1)	1	0
B(1)	2	0
B(1)	3	0
B(6)    2   0
B(6)    3   0
B(8)    3   0
];
%% Assemblage de la force exterieure
% global DISPTD
% [NUMNP, NDOF] = size(XYZ);				
% NE = size(LE,1);
% NEQ = NDOF*NUMNP;
% FORCEXTG = sparse(NEQ+NE,1); % FORCE EXTERIEUR GENERALISEE
% DISPTD= sparse(NEQ+NE,1);
% if size(EXTFORCE,1)>0
% LOC = NDOF*(EXTFORCE(:,1)-1)+EXTFORCE(:,2);
% FORCEXTG(LOC) = FORCEXTG(LOC) + EXTFORCE(:,3);
% end

% Set program parameters
ITRA=100; ATOL=1.0E5; NTOL=6; TOL=1E-6;
% Load increments [Start End Increment InitialFactor FinalFactor]
TIMS=[0.0 1 1 0.0 1]';
MID=0;
% TOL=1.0e-6;
% itermax=100; 
% iter=0;
% 
% solveur(TOL,itermax,FORCEXTG,SDISPT,PROP, XYZ,LE)
NOUT = fopen('flexioncoin.txt','w');
lineSearch=false;
NLFEA2(ITRA,TOL,ATOL,NTOL,TIMS,NOUT,MID,PROP,EXTFORCE,SDISPT,XYZ,LE,lineSearch)