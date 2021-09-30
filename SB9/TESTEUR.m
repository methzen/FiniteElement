clc
clear all

%%
 XYZ=[0 0 0
      2 0 0;
      2 2 0;
      0 2 0;
      0 0 2;
      2 0 2;
      2 2 2;
      0 2 2];
LE=[1 2 3 4 5 6 7 8];
E=1;
v=0.3;
PROP=[E   v];
%% Conditions limites
% External forces [Node, DOF, Value]
F=-1;
A=LE(end,[1 2 6 5]);
%A=LE(end,[2 3 6 7]);
%A=LE(end,[5 6 7 8]);
EXTFORCE=[A(1) 2 F
          A(2) 2 F
          A(3) 2 F
          A(4) 2 F];

% Prescribed displacements [Node, DOF, Value]

B=LE(1,[3 4 7 8]);
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
];

%%
ITRA=1000; 
ATOL=1.0E10; 
NTOL=6; 
TOL=1E-6;

% Load increments [Start End Increment InitialFactor FinalFactor]
TIMS=[0.0 1 0.5 0.0 1]';
MID=0;
%
NOUT = fopen('cube.txt','w');
lineSearch=false;
NLFEA2(ITRA,TOL,ATOL,NTOL,TIMS,NOUT,MID,PROP,EXTFORCE,SDISPT,XYZ,LE,lineSearch)

















% [BG0,BGxi,BGeta,BGzeta,BGxieta,BGxizeta,BGetazeta]=matriceBg(XYZ);
% [B0,BZETA,BZZ,BXI,BETA,BETAZETA,BXIZETA]=matriceBcart(XYZ,xyz);
% E3=0;


% E0=deformEcart(DSP,XYZ,[0 0 E3])
% B9=calculB9(XYZ,[0 0 E3]);
% 
%       
% BG=BG0+E3*BGzeta;
% Bm=B0+E3*BZETA+E3*E3*BZZ;
% [T0,TXI,TETA,TZETA]=matriceTdecomp(XYZ) ;
% 

% xyz=XYZ+D';
% T=matriceT(XYZ,[0 0 0])
% [Ec0,Eczeta,Eczz,EcXI,EcETA,EcETAZETA,EcXIZETA]=deformatEc(D,XYZ)
% 
% [Bc0,BcZETA,BcZZ,BcXI,BcETA,BcETAZETA,BcXIZETA]=matriceBc(XYZ)
% [B0,BZETA,BZZ,BXI,BETA,BETAZETA,BXIZETA]=matriceBcart(XYZ,xyz)
% XI=[0 0 0];
%  STRESS=[1 1 1 0 0 0]';
%  Kg=rigigeo(STRESS,XYZ,XI);
%  [EKF,K]=Kgeo3d(XYZ, STRESS, XI);
%  EKF;
%  
%           % calcul de la rigidité géo
%          SIG=[STRESS(1) STRESS(4) STRESS(6);
%               STRESS(4) STRESS(2) STRESS(5);
%               STRESS(6) STRESS(5) STRESS(3)];
%         SHEAD=zeros(9);
%         SHEAD(1:3,1:3)=SIG;
%         SHEAD(4:6,4:6)=SIG;
%         SHEAD(7:9,7:9)=SIG;
% [BG0,BGxi,BGeta,BGzeta,BGxieta,BGxizeta,BGetazeta]=matriceBg(XYZ);
% BG=BG0+E3*BGzeta; 
% Km=BG'*SHEAD*BG;
% [B0,BZETA,BZZ,BXIdev,BETAdev,BETAZETAdev,BXIZETAdev]=matriceBcart(XYZ,xyz)