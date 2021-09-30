function [Bc0,Bczeta,Bczz,BcXI,BcETA,BczetaXI,BczetaETA,BcXIXI,BcETET,BcXIETA]=matriceBc(XYZ)

% ICI LES COORDONNÉES SONT LES COORDONNÉES INITIALES
%  ALERTE ECRITURE VERIFIEE
%
[J10,J1ZETA,...
 J20,J2ZETA,...
 J30,J3ETA,J3XI]=jacobienDecompose(XYZ);

% initialisation :
[g1,g2,g3,h1,h2]=g1g2g3h1h2;
BczetaXI=zeros(6,18);
BczetaETA=zeros(6,18);
BcXIETA=zeros(6,18);
BcETET=zeros(6,18);
BcXIXI=zeros(6,18);
BcETA=zeros(6,18);
BcXI=zeros(6,18);
Bczz=zeros(6,18);
Bczeta=zeros(6,18);
Bc0=zeros(6,18);
%% TYING POINTS FOR TS
for I=1:6
COL=3*(I-1)+1:3*(I-1)+3;
XIA1=[0.5 0 -1];
XIA2=[0.5 0 +1];
XIB1=[0 0.5 -1];
XIB2=[0 0.5 +1];
XIC1=[0.5 0.5 -1];
XIC2=[0.5 0.5 +1];
%
[E13A1,E23A1]=Bct(XYZ,XIA1,I);
[E13B1,E23B1]=Bct(XYZ,XIB1,I);
[E13C1,E23C1]=Bct(XYZ,XIC1,I);
%
[E13A2,E23A2]=Bct(XYZ,XIA2,I);
[E13B2,E23B2]=Bct(XYZ,XIB2,I);
[E13C2,E23C2]=Bct(XYZ,XIC2,I);
%
% FACE INF
E23_01=E23B1;
E23_xi1=-(E23B1-E13A1-E23C1+E13C1);
%
E13_01=E13A1;
E13_eta1=(E23B1-E13A1-E23C1+E13C1);
% FACE SUP
E23_02=E23B2;
E23_xi2=-(E23B2-E13A2-E23C2+E13C2);
%
E13_02=E13A2;
E13_eta2=(E23B2-E13A2-E23C2+E13C2);

%%
Bc0(1,COL)=J10*g1(I);
Bc0(2,COL)=J20*g2(I);
Bc0(3,COL)=J30*g3(I);
Bc0(4,COL)=J10*g2(I)+J20*g1(I);
Bc0(5,COL)=0.5*(E13_02+E13_01);
Bc0(6,COL)=0.5*(E23_02+E23_01);

%%
Bczeta(1,COL)=J10*h2(I)+g1(I)*J1ZETA;
Bczeta(2,COL)=J20*h1(I)+g2(I)*J2ZETA;
Bczeta(3,COL)=[0 0 0];
Bczeta(4,COL)=J10*h1(I)+J1ZETA*g2(I)+J20*h2(I)+J2ZETA*g1(I);
Bczeta(5,COL)=0.5*(E13_02-E13_01);
Bczeta(6,COL)=0.5*(E23_02-E23_01);
%%
Bczz(1,COL)=J1ZETA*h2(I);
Bczz(2,COL)=J2ZETA*h1(I);
Bczz(4,COL)=J1ZETA*h1(I)+J2ZETA*h2(I);
%% 
BcXI(3,COL)=J30*h2(I)+J3XI*g3(I);
BcXI(5,COL)=[0 0 0];
BcXI(6,COL)=0.5*(E23_xi2+E23_xi1);
%%
BcETA(3,COL)=J30*h1(I)+J3ETA*g3(I);
BcETA(5,COL)=0.5*(E13_eta1+E13_eta2);
BcETA(6,COL)=0;
%%
BcXIXI(3,COL)=J3XI*h2(I);
%%
BcETET(3,COL)=J3ETA*h1(I);
%%
BcXIETA(3,COL)=(J3XI*h1(I))+(h2(I)*J3ETA);
%%
BczetaETA(5,COL)=0.5*(E13_eta2-E13_eta1);
%%
BczetaXI(6,COL)=0.5*(E23_xi2-E23_xi1);
end

%% Changement de base
Bc0=Bc0-(BcXI+BcETA)/3+(BcXIXI+BcETET+BcXIETA)/9;
Bczeta=Bczeta-(BczetaXI+BczetaETA)/3;
BcXI=BcXI-2*BcXIXI/3-BcXIETA/3;
BcETA=BcETA-2*BcETET/3-BcXIETA/3;
end