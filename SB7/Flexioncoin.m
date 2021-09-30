%
% Two-element example page 131
%
close all
clear all
clc
format long g
%% Donnee d'entrees
% Nodal coordinates
% geom=Plaque(100,100,1,1,1,1);
% XYZ=full(geom.Nodes3D(:,2:4));
% % %
% % % Element connectivity
%LE=geom.Connect3D;
%LE=[1 2 3 4 5 6 7 8];
%
 XYZ=[ 0 0 0;
       100 0 0;
       100 100 0;
       0 100 0;
       0 0 1;
       100 0 1;
       100 100 1;
       0 100 1];
LE=[1 2 3 4 5 6 7 8];
%LE=Mailcubeenprisme(LE);
[LE,XYZ]=Mailcubeenprisme(LE,XYZ,1);

%  LE=[2 4 1 6 8 5;
%      4 2 3 8 6 7];
% %   B=gradtransf(XYZ,[0 0 0]),, 
for i=1:size(LE,1)
plotHexa6(XYZ(LE(i,:),:),'b')
end
%AfficheRepere(0,0,0,eye(3),0)
%% Material properties HYPERELASTIQUE
%
E=200000;
v=0.3;
PROP=[E   v];
%% Conditions limites
% External forces [Node, DOF, Value]
F=1;
%B=LE(1,[1 2 3 4 5 6 7 8]);
EXTFORCE=[3 3 F];
% Prescribed displacements [Node, DOF, Value]
SDISPT=[
1	1	0
1	2	0
1	3	0
6   2   0
6   3   0
8   3   0];
%% Assemblage de la force exterieure
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