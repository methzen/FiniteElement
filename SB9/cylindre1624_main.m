%
% Two-element example page 131
%
clear all
clc
format long g
%% Donnee d'entrees
% Nodal coordinates
[XYZ,LE,cx,cy,cz]=meshCylind1624;
%geom=Plaque(10,1,0.1,16,1,1);
%XYZ=full(geom.Nodes3D(:,2:4));
% % % Element connectivity
%LE=geom.Connect3D;
for i=1:size(LE,1)
plotHexa8(XYZ(LE(i,:),:),'b')
end
%AfficheRepere(0,0,0,eye(3),0)
%% Material properties HYPERELASTIQUE
%
E=10.5e6;
v=0.3125;
PROP=[E   v];
%% Conditions limites
% External forces [Node, DOF, Value]
F=10000;
EXTFORCE=[58 2 F
          850 2 F];
% Prescribed displacements [Node, DOF, Value]

SDISPT=[cx;cy;cz];

%%
ITRA=1000; 
ATOL=1.0E10; 
NTOL=6; 
TOL=1E-6;

% Load increments [Start End Increment InitialFactor FinalFactor]
TIMS=[0.0 1 0.01 0.0 1]';
MID=-1;
%
NOUT = fopen('cylindre1624.txt','w');
lineSearch=false;
NLFEA2(ITRA,TOL,ATOL,NTOL,TIMS,NOUT,MID,PROP,EXTFORCE,SDISPT,XYZ,LE,lineSearch)