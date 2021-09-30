
clear all
clc
format long g
%% Geometrie importation un maillage plan
[XYZ,LE,g1,g2,g3]=Cylindre_Salome; 

%% Deplacement impose
SDISPT=[
 g1    ones(size(g1)) zeros(size(g1));
 g2  2*ones(size(g2)) zeros(size(g2));
 g3  3*ones(size(g3)) zeros(size(g3))];
%% Chargement
% External forces [Node, DOF, Value]
F=-50;
EXTFORCE=[1 1 F
          2 1 F];

%% materiau
E=10.5e6;
v=0.3125;
PROP=[E   v];
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
% %
% TOL=1.0e-6;
% itermax=100; 
% %
% solveur(TOL,itermax,FORCEXTG,SDISPT,PROP, XYZ,LE)
% 
% full(DISPTD(1:6))

%% d�form�e calcul Elem Fini


%ref=0.1139
%bassa=0.11206

ITRA=70; ATOL=1.0E5; NTOL=6; TOL=1E-6;
% Load increments [Start End Increment InitialFactor FinalFactor]
TIMS=[0.0 1 1 0.0 1]';
MID=0;
% TOL=1.0e-6;
% itermax=100; 
% iter=0;
% 
% solveur(TOL,itermax,FORCEXTG,SDISPT,PROP, XYZ,LE)
NOUT = fopen('cylindrepincee.txt','w');
lineSearch=false;
NLFEA2(ITRA,TOL,ATOL,NTOL,TIMS,NOUT,MID,PROP,EXTFORCE,SDISPT,XYZ,LE,lineSearch)

