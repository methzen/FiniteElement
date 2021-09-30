%
% Two-element example page 131
%
clear all
clc
format long g
%% Donnee d'entrees
% Nodal coordinates
geom=Plaque(10,1,0.1,16,1,1);
XYZ=full(geom.Nodes3D(:,2:4));
% % % Element connectivity
LE=geom.Connect3D;
for i=1:size(LE,1)
plotHexa8(XYZ(LE(i,:),:),'b')
end
%AfficheRepere(0,0,0,eye(3),0)
%% Material properties HYPERELASTIQUE
%
E=1e9;
v=0.3;
PROP=[E   v];
%% Conditions limites
% External forces [Node, DOF, Value]
F=50;
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
B(2)   1	0
B(2)   2   0
B(2)   3   0
B(3)   1	0
B(3)   2   0
B(3)   3   0
B(4)	1	0
B(4)	2	0
B(4)	3	0
];

%% Assemblage de la force exterieure
global DISPTD
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
%%
CON.nincr          = 0;
CON.xlmax          = fscanf(fid,'%g',1);
CON.dlamb          = fscanf(fid,'%g',1);
CON.miter          = fscanf(fid,'%d',1);
CON.cnorm          = fscanf(fid,'%g',1);
CON.searc          = fscanf(fid,'%g',1);
CON.xlamb = 0;
CON.incrm = 0; 
CON.searc        = 1e5; 
CON.msearch         = 5;
CON.incrm           = 0;
LSNR(CON,FORCEXTG,SDISPT,PROP,XYZ,LE)
%