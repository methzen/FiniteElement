%
% Two-element example page 131
%
clear all
clc
format long g
%% Donnee d'entrees
% Nodal coordinates
geom=Plaque(1,1,0.1,4,4,1);
XYZ=full(geom.Nodes3D(:,2:4));
% % % Element connectivity
LE=geom.Connect3D;
%LE=Mailcubeenprisme(LE);
[LE,XYZ]=Mailcubeenprisme(LE,XYZ,1);
for i=1:size(LE,1)
plotHexa6(XYZ(LE(i,:),:),'b')
end
%AfficheRepere(0,0,0,eye(3),0)
%% Material properties HYPERELASTIQUE
%
E=210e9;
v=0.;
PROP=[E   v];
%% Conditions limites
% External forces [Node, DOF, Value]
F=100;
%A=LE(end,[2 3 6 7]);
A=[21 22 23 24 25 46 47 48 49 50];
EXTFORCE=[A(1) 3 F
          A(2) 3 F
          A(3) 3 F
          A(4) 3 F];

% Prescribed displacements [Node, DOF, Value]

%B=LE(1,[1 4 5 8]);
SDISPT=zeros(30,3);
B=[1 2 3 4 5 26 27 28 29 30];
for i=1:10
    lin=3*(i-1)+1:3*(i-1)+3;
SDISPT(lin,:)=[B(i)	1	0;
B(i)	2	0;
B(i)	3	0];
end

%%
ITRA=3; 
ATOL=1.0E10; 
NTOL=1; 
TOL=1E-5;
% Load increments [Start End Increment InitialFactor FinalFactor]
TIMS=[0.0 1 1 0.0 1]';
MID=0;
%
NOUT = fopen('poutretest.txt','w');
lineSearch=false;
NLFEA2(ITRA,TOL,ATOL,NTOL,TIMS,NOUT,MID,PROP,EXTFORCE,SDISPT,XYZ,LE,lineSearch)