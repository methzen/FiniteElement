
function [T0,TXI,TETA,TZETA]=matriceTdecomp(XYZ)
%
X=XYZ(:,1);
Y=XYZ(:,2);
Z=XYZ(:,3);
%
[g1,g2,g3,h1,h2]=g1g2g3h1h2;
%
J0=Jacobien(XYZ,[1/3 1/3 0]);
J0=inv(J0);
%
JXI=[0 0 dot(X,h2);
     0 0 dot(Y,h2);
     0 0 dot(Z,h2)];
%    
JETA=[0 0 dot(X,h1);
      0 0 dot(Y,h1);
      0 0 dot(Z,h1)];
%  
JZETA=[dot(X,h2) dot(X,h1) 0;
       dot(Y,h2) dot(Y,h1) 0;
       dot(Z,h2) dot(Z,h1) 0];
%
A=-J0*JXI*J0;
B=-J0*JETA*J0;
C=-J0*JZETA*J0;  
%%
T0=matriceT(XYZ,[1/3 1/3 0]);
[TXI,TETA,TZETA]=matricetxi(J0,A,B,C);

end
   