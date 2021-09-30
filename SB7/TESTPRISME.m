clear all
close all
clc  
%
 XYZ=[0     0    -1
      1     0    -1
      0     1    -1
      0     0     1
      1     0     1
      0     1     1
      ];
LE=[1 2 3 4 5 6];
%Nu=zeros(1,18);
%PrintMesh(Nu,XYZ, Ec,1)
[D,C]= matrice_D(1,0.3);
%
ZG=[-1 -sqrt(3/7) 0 sqrt(3/7)  1];
WGT=[0.1 49.0/90.0 32.0/45.0 49.0/90.0 0.1];
%% 
ni=1;
nj=2;
nk=3;
nddl=18;
Be=zeros(6,nddl);
matDint=zeros(3);
[DTAN,~]= matrice_D(1, 0.3);
Nijk=[1 2 3];
 K=zeros(nddl,nddl);
for i=1:5
J=Jacobien(XYZ,[1/3 1/3 ZG(i)]);
FAC=det(J)*WGT(i);        
Be=calculBcartassume(XYZ,XYZ,[1/3 1/3 ZG(i)]);
K=K+Be'*DTAN*Be*FAC;
end
Kmx=K;
%%
 Nu=null(Kmx);
 size(Nu)
%%
DSP1=[
      1 -1  0
      1  1  0
     -1 -1  0
      1 -1  0
      1  1  0
     -1 -1  0
]';

DSP2 =[
       2 0 -1 
       2 0  1  
       2 0 -1
      -2 0 -1
      -2 0  1  
      -2 0 -1
]';

DSP3 =[0  2 -1 
       0  2 -1  
       0  2  1
       0 -2 -1      
       0 -2 -1  
       0 -2  1
 ]';
  
DSP4 =[0  0  1
       0  0  1  
       0  0  1
       0  0  1 
       0  0  1  
       0  0  1]';
   
DSP5 =[0  1  0
       0  1  0  
       0  1  0
       0  1  0 
       0  1  0 
       0  1  0]';
   
DSP6 =[1  0  0
       1  0  0  
       1  0  0
       1  0  0 
       1  0  0  
       1  0  0]';
   
DSP7 =[0  0  -1
       0  0  +2  
       0  0  -1
       0  0  +1 
       0  0  -2  
       0  0  +1]';
   
DSP8 =[0  0  -1
       0  0  -1  
       0  0  +2
       0  0  +1 
       0  0  +1  
       0  0  -2]'; 
   
T=[ -1
     1
     0
    -1
    -2
     0
     2
     1
     0
     1
    -1
     0
     1
     2
     0
    -2
    -1
     0];
V=Nu;
DSP9=V(:,1)+V(:,2)+V(:,3)+V(:,4)+V(:,5)+V(:,6)+V(:,7);
A=[DSP4(:) DSP5(:) DSP6(:) DSP1(:) DSP2(:) DSP3(:) DSP9];
Q=gramshmid(A);

for i=1:7
    PrintMesh(Q(:,i),XYZ, LE,1); 
end