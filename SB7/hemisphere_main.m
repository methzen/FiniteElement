%
% Two-element example page 131
%
clear all
clc
format long g
%% Donnee d'entrees
% Nodal coordinates
[XYZ,LE,Cx,Cy]=hemisphere;
[LE,XYZ]=Mailcubeenprisme(LE,XYZ,1);
%geom=Plaque(10,1,0.1,16,1,1);
%XYZ=full(geom.Nodes3D(:,2:4));
% % % Element connectivity
%LE=geom.Connect3D;
for i=1:size(LE,1)
plotHexa6(XYZ(LE(i,:),:),'b')
end
%AfficheRepere(0,0,0,eye(3),0)
%% Material properties HYPERELASTIQUE
%
E=6.825e7;
v=0.3;
PROP=[E   v];
%% Conditions limites
% External forces [Node, DOF, Value]
F=0.5;
EXTFORCE=[1 1 F
          2 1 F
          50 2 -F
          66 2 -F];
% Prescribed displacements [Node, DOF, Value]

SDISPT=[Cx;Cy;1 3 0];

%%
ITRA=1000; 
ATOL=1.0E10; 
NTOL=6; 
TOL=1E-5;

% Load increments [Start End Increment InitialFactor FinalFactor]
TIMS=[0.0 1 1 0.0 1]';
MID=0;
%
NOUT = fopen('hemisphere.txt','w');
lineSearch=false;
NLFEA2(ITRA,TOL,ATOL,NTOL,TIMS,NOUT,MID,PROP,EXTFORCE,SDISPT,XYZ,LE,lineSearch)