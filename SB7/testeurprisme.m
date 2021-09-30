
close all
clc
% XYZ=[0     0    -0.1
%      2     0    -0.1
%      2     1    -0.1
%      0     1    -0.1
%      1   0.5    -0.1
%      0     0    +0.1
%      2     0    +0.1
%      2     1    +0.1
%      0     1    +0.1
%      1   0.5    +0.1];

XYZ=[1     0    -1
     0     1    -1
     0     0    -1    
     1     0    +1
     0     1    +1
     0     0    +1];
%
LE=[1 2 3 4 5 6];



% T=matriceT(XYZ,[0 0 0]);
% E0=T*Ec0
% DSP(:);
% B0=T*Bc0*DSP(:);

% LE=[1 2 5 6 7 10
%     2 3 5 7 8 10
%     3 4 5 8 9 10
%     4 1 5 9 6 10];

% XYZ=XYZ(LE(1,:),:);
% 
%  XI=[1/6 1/6 -1/sqrt(3);
%      1/6 2/3 -1/sqrt(3);
%      2/3 1/6 -1/sqrt(3);
%      1/6 1/6 +1/sqrt(3);
%      1/6 2/3 +1/sqrt(3);
%      2/3 1/6 +1/sqrt(3)];
% Benh=zeros(3); 
nddl=18;
 K=zeros(nddl);
  % Integration points and weights
ZG=[-1 -sqrt(3/7) 0 sqrt(3/7)  1];
WGT=[0.1 49.0/90.0 32.0/45.0 49.0/90.0 0.1];


%             ZG(1) = -0.906179845938664d0;
%             ZG(2) = -0.538469310105683d0;
%             ZG(3) = 0.d0;
%             ZG(4) = 0.538469310105683d0;
%             ZG(5) = 0.906179845938664d0;
% 
%             WGT(1) = 0.236926885056189d0;
%             WGT(2) = 0.478628670499366d0;
%             WGT(3) = 0.568888888888889d0;
%             WGT(4) = 0.478628670499366d0;
%             WGT(5) = 0.236926885056189d0;


g1=[1 0 -1 1 0 -1]'/2;
g2=[0 1 -1 0 1 -1]'/2;
g3=[0 0 -1 0 0  1]'/2;
h1=[-1 0 1 1 0 -1]'/2;
h2=[0 -1 1 0 1 -1]'/2;
%N=[r*(1-t) s*(1-t) (1-r-s)*(1-t) r*(1+t) s*(1+t) (1-r-s)*(1+t)];



BN=zeros(6,nddl);
%     % Determinant and shape function derivatives
for LZ=1:size(WGT,2);
      %E3=ZG(LZ);
      J0=Jacobien(XYZ,[1/3 1/3 ZG(LZ)]);
      DET= det(J0);
      FAC=WGT(LZ)*DET/2;
      BN=calculBcartassume(XYZ,[1/3 1/3 ZG(LZ)]);      
      [DTAN,~]= matrice_D(10, 0.3);
      % Loi de comportement de Saint Venant Kirchhoff
      % Tangent stiffness
      K=K+BN'*DTAN*BN*FAC;
end
%     
size(null(K),2)

DSP1=[
      1  1  0
     -1 -1  0
      1 -1  0
      1  1  0
     -1 -1  0
      1 -1  0]';

DSP2 =[2 0  1  
       2 0 -1
       2 0 -1 
      -2 0  1  
      -2 0 -1
      -2 0 -1]';

DSP3 =[0  2 -1  
       0  2  1
       0  2 -1       
       0 -2 -1  
       0 -2  1
       0 -2 -1 ]';
  
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
% K*DSP(:);
%PrintMesh(DSP1(:),XYZ, LE,1);

V=null(K);
DSP7=V(:,1)+V(:,2)+V(:,3)+V(:,4)+V(:,5)+V(:,6);
A=[DSP4(:) DSP5(:) DSP6(:) DSP1(:) DSP2(:) DSP3(:) DSP7];
Q=gramshmid(A);

for i=1:7
    U=V(:,i);
    M=U'*K*U;
    %PrintMesh(Q(:,i),XYZ, LE,1); 
end
DSPZ =[ 0  1 -1
       -1  0 -1  
        0  0 -1
        0  1  1 
       -1  0  1  
        0  0  1]';
    
T=[ -1
    -2
     0
     2
     1
     0
    -1
     1
     0
     1
     2
     0
    -2
    -1
     0
     1
    -1
     0];
[Bc,Br,Bs,Bz,Bzz,Brz,Bsz,Brs]=matriceBcart(XYZ,XYZ);
 