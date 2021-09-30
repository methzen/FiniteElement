%
% Two-element example page 131
%
clear all
clc
format long g
%% Donnee d'entrees
% Nodal coordinates
[XYZ,LE]=annularplatedata ;
% % % Element connectivity
for i=1:size(LE,1)
plotHexa8(XYZ(LE(i,:),:),'b')
end
% %AfficheRepere(0,0,0,eye(3),0)
for i=1:size(LE,1)
plotHexa8(XYZ(LE(i,:),:),'r')
end
%% Material properties HYPERELASTIQUE
%
E=21e6;
v=0.0;
PROP=[E   v];
%% Conditions limites
% External forces [Node, DOF, Value]
alpha=0.8*2/24;
A=[1 2 3 4 5 6 7 8 9 10 11 12 13 14 ]';
dir=2;
extr=ones(1,4)*alpha;
inter=ones(1,10)*2*alpha;
F=[extr inter]';
EXTFORCE=[A dir*ones(size(A)) F];

% Prescribed displacements [Node, DOF, Value]

B=[74 44 104 134 164 194 224 254 284 314 344 374 404 434];
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
B(4)	3	0
B(5)	1	0
B(5)	2	0
B(5)	3	0
B(6)	1	0
B(6)	2	0
B(6)	3	0
B(7)	1	0
B(7)	2	0
B(7)	3	0
B(8)    1	0
B(8)    2   0
B(8)    3   0
B(9)	1	0
B(9)    2   0
B(9)    3   0
B(10)	1	0
B(10)	2	0
B(10)	3	0
B(11)	1	0
B(11)	2	0
B(11)	3	0
B(12)	1	0
B(12)	2	0
B(12)	3	0
B(13)	1	0
B(13)	2	0
B(13)	3	0
B(14)	1	0
B(14)	2	0
B(14)	3	0
];

%%
ITRA=1000; 
ATOL=1.0E9; 
NTOL=6; 
TOL=1E-6;

% Load increments [Start End Increment InitialFactor FinalFactor]
TIMS=[0.0 1 0.05 0.0 1]';
MID=-1;
%
NOUT = fopen('annularplate.txt','w');
lineSearch=false;
NLFEA2(ITRA,TOL,ATOL,NTOL,TIMS,NOUT,MID,PROP,EXTFORCE,SDISPT,XYZ,LE,lineSearch)