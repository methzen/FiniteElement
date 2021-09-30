
function [BG0,BGxi,BGeta,BGzeta,BGxieta,BGxizeta,BGetazeta]=matriceBg(XYZ)

%%
X=XYZ(:,1);
Y=XYZ(:,2);
Z=XYZ(:,3);
%
g1=[-1 1 1 -1 -1 1 1 -1]/8;
g2=[-1 -1 1 1 -1 -1 1 1]/8;
g3=[-1 -1 -1 -1 1 1 1 1]/8;
%
h1=[1 -1 1 -1 1 -1 1 -1]/8;
h2=[1 1 -1 -1 -1 -1 1 1]/8;
h3=[1 -1 -1 1 -1 1 1 -1]/8;
h4=[-1 1 -1 1 1 -1 1 -1]/8;
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
%
A=A';
B=B';
C=C';
J0=J0';


Nc0=[g1;g2;g3];
Ncxi=[zeros(1,8);h1;h3];
Nceta=[h1;zeros(1,8);h2];
Nczeta=[h3;h2;zeros(1,8)];
Ncxieta=[zeros(1,8);zeros(1,8);h4];
Ncxizeta=[zeros(1,8);h4;zeros(1,8)];
Ncetazeta=[h4;zeros(1,8);zeros(1,8)];
%
N0=J0*Nc0;
Nxi=A*Nc0+J0*Ncxi;
Neta=B*Nc0+J0*Nceta;
Nzeta=C*Nc0+J0*Nczeta;
Nxieta=J0*Ncxieta;
Nxizeta=J0*Ncxizeta;
Netazeta=J0*Ncetazeta;
%
BG0=zeros(9,24);
BGxi=zeros(9,24);
BGeta=zeros(9,24);
BGzeta=zeros(9,24);
BGZZ=zeros(9,24);
BGxieta=zeros(9,24);
BGxizeta=zeros(9,24);
BGetazeta=zeros(9,24);
for I=1:8
    COL=(I-1)*3+1:(I-1)*3+3;
       BG0(:,COL)=[N0(1,I) 0         0;
                   N0(2,I) 0         0;
                   N0(3,I) 0         0;
                   0         N0(1,I) 0;
                   0         N0(2,I) 0;
                   0         N0(3,I) 0;
                   0         0         N0(1,I);
                   0         0         N0(2,I);
                   0         0         N0(3,I)];
%               
       BGxi(:,COL)=[Nxi(1,I) 0         0;
                    Nxi(2,I) 0         0;
                    Nxi(3,I) 0         0;
                    0         Nxi(1,I) 0;
                    0         Nxi(2,I) 0;
                    0         Nxi(3,I) 0;
                    0         0         Nxi(1,I);
                    0         0         Nxi(2,I);
                    0         0         Nxi(3,I)]; 
%                
      BGeta(:,COL)=[Neta(1,I) 0         0;
                    Neta(2,I) 0         0;
                    Neta(3,I) 0         0;
                    0         Neta(1,I) 0;
                    0         Neta(2,I) 0;
                    0         Neta(3,I) 0;
                    0         0         Neta(1,I);
                    0         0         Neta(2,I);
                    0         0         Neta(3,I)];
%                
     BGzeta(:,COL)=[Nzeta(1,I) 0         0;
                    Nzeta(2,I) 0         0;
                    Nzeta(3,I) 0         0;
                    0         Nzeta(1,I) 0;
                    0         Nzeta(2,I) 0;
                    0         Nzeta(3,I) 0;
                    0         0         Nzeta(1,I);
                    0         0         Nzeta(2,I);
                    0         0         Nzeta(3,I)]; 
%                
    BGxieta(:,COL)=[Nxieta(1,I) 0         0;
                    Nxieta(2,I) 0         0;
                    Nxieta(3,I) 0         0;
                    0         Nxieta(1,I) 0;
                    0         Nxieta(2,I) 0;
                    0         Nxieta(3,I) 0;
                    0         0         Nxieta(1,I);
                    0         0         Nxieta(2,I);
                    0         0         Nxieta(3,I)];
%                
   BGxizeta(:,COL)=[Nxizeta(1,I) 0         0;
                    Nxizeta(2,I) 0         0;
                    Nxizeta(3,I) 0         0;
                    0         Nxizeta(1,I) 0;
                    0         Nxizeta(2,I) 0;
                    0         Nxizeta(3,I) 0;
                    0         0         Nxizeta(1,I);
                    0         0         Nxizeta(2,I);
                    0         0         Nxizeta(3,I)];
 %               
  BGetazeta(:,COL)=[Netazeta(1,I) 0         0;
                    Netazeta(2,I) 0         0;
                    Netazeta(3,I) 0         0;
                    0         Netazeta(1,I) 0;
                    0         Netazeta(2,I) 0;
                    0         Netazeta(3,I) 0;
                    0         0         Netazeta(1,I);
                    0         0         Netazeta(2,I);
                    0         0         Netazeta(3,I)];               
end
end