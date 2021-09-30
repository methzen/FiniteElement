%
% Two-element example page 131
%
clear all
clc
format long g
%% Donnee d'entrees
% Nodal coordinates
[XYZ,LE]=tubedistgeom();
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
E=1e4;
v=0.3;
PROP=[E   v];
%% Conditions limites
% External forces [Node, DOF, Value]
F=10;
A=LE(end,[2 3 6 7]);
EXTFORCE=[A(1) 3 F
          A(2) 3 F
          A(3) 3 F
          A(4) 3 F];

% Prescribed displacements [Node, DOF, Value]


B=LE(1,[1 4 5 8]);
SDISPT=[B(1)	1	0
B(1)	2	0
B(1)	3	0
B(2)    1	0
B(2)    2   0
B(2)    3   0
B(3)	1	0
B(3)    2   0
B(3)    3   0
B(4)	1	0
B(4)	2	0
B(4)	3	0];

%%
ITRA=1; 
ATOL=1.0E10; 
NTOL=1; 
TOL=1E-6;

% Load increments [Start End Increment InitialFactor FinalFactor]
TIMS=[0.0 1 1 0.0 1]';
MID=-1;
%
NOUT = fopen('tubedis.txt','w');
lineSearch=false;
NLFEA2(ITRA,TOL,ATOL,NTOL,TIMS,NOUT,MID,PROP,EXTFORCE,SDISPT,XYZ,LE,lineSearch)