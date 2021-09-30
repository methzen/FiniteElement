
function [T0,TXI,TETA,TZETA]=matriceTdecomp(XYZ)
%
X=XYZ(:,1);
Y=XYZ(:,2);
Z=XYZ(:,3);
%
h1=[1 -1 1 -1 1 -1 1 -1]/8;
h2=[1 1 -1 -1 -1 -1 1 1]/8;
h3=[1 -1 -1 1 -1 1 1 -1]/8;
%
J0=Jacobien(XYZ,[0 0 0]);
J0=inv(J0);
%
JXI=[0 dot(X,h1) dot(X,h3);
     0 dot(Y,h1) dot(Y,h3);
     0 dot(Z,h1) dot(Z,h3)];
%    
JETA=[dot(X,h1) 0 dot(X,h2);
      dot(Y,h1) 0 dot(Y,h2);
      dot(Z,h1) 0 dot(Z,h2)];
%  
JZETA=[dot(X,h3) dot(X,h2) 0;
       dot(Y,h3) dot(Y,h2) 0;
       dot(Z,h3) dot(Z,h2) 0];
%
A=-J0*JXI*J0;
B=-J0*JETA*J0;
C=-J0*JZETA*J0;  
%%
T0=matriceT(XYZ,[0 0 0]);
[TXI,TETA,TZETA]=matricetxi(J0,A,B,C);
%  TXI=matriceTest(A);
%  TETA=matriceTest(B);
%  TZETA=matriceTest(C);

end
   