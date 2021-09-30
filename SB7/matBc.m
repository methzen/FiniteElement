function Bc0=matBc(XYZ,XI)

% ICI IL FAUT RENTRER AVEC LES POSITIONS ACTUELLES 
[J1,J2,J3]=colomnejacobien(XYZ,XI);
      
%%
[g1,g2,g3,h1,h2]=g1g2g3h1h2;
%
%% TYING POINTS FOR CT AND PINCH
d=1/100000;
XIA1=[1/6 1/3 -1];
XIA2=[1/6 1/3 +1];
XIB1=[1/3 1/6 -1];
XIB2=[1/3 1/6 +1];
XIC1=[1/6 1/6 -1];
XIC2=[1/6 1/6 +1];
XID1=[4*d d -1];
XID2=[4*d d +1];
XIE1=[d 4*d -1];
XIE2=[d 4*d +1];
XIF1=[4*d 4*d -1];
XIF2=[4*d 4*d +1];
%
XIA=[1 0 0];
XIB=[0 1 0];
XIC=[0 0 0];
%
N1=g1+XI(3)*h1;
N2=g2+XI(3)*h2;
N3=g3+XI(1)*h1+XI(2)*h2;
%%
 Bc0=zeros(6,18);
for I=1:6
%%
[B13A1,B23A1]=Bct(XYZ,XIA1,I);
[B13B1,B23B1]=Bct(XYZ,XIB1,I);
[B13C1,B23C1]=Bct(XYZ,XIC1,I);
[B13D1,B23D1]=Bct(XYZ,XID1,I);
[B13E1,B23E1]=Bct(XYZ,XIE1,I);
[B13F1,B23F1]=Bct(XYZ,XIF1,I);
[B13A2,B23A2]=Bct(XYZ,XIA2,I);
[B13B2,B23B2]=Bct(XYZ,XIB2,I);
[B13C2,B23C2]=Bct(XYZ,XIC2,I);
[B13D2,B23D2]=Bct(XYZ,XID2,I);
[B13E2,B23E2]=Bct(XYZ,XIE2,I);
[B13F2,B23F2]=Bct(XYZ,XIF2,I);
B33A=Bpincementcovari(XYZ,XIA,I);
B33B=Bpincementcovari(XYZ,XIB,I);
B33C=Bpincementcovari(XYZ,XIC,I);
%
%FACE INF
b23_01=(2*B23A1-B13A1+B13C1+B23C1+B13F1-B23F1-B13D1+B23E1)/3.0;
b23_xi1=B13D1+B23F1-B13F1-B23E1;
b13_01=(2*B13B1-B23B1+B13C1+B23C1-B13F1+B23F1+B13D1-B23E1)/3.0;
b13_eta1=B13F1-B13D1-B23F1+B23E1;
%Facesup
b23_02=(2*B23A2-B13A2+B13C2+B23C2+B13F2-B23F2-B13D2+B23E2)/3.0;
b23_xi2=B13D2+B23F2-B13F2-B23E2;
b13_02=(2*B13B2-B23B2+B13C2+B23C2-B13F2+B23F2+B13D2-B23E2)/3.0;
b13_eta2=B13F2-B13D2-B23F2+B23E2;
%
b23=(b23_02-b23_01+XI(1)*(b23_xi2-b23_xi1))*XI(3)/2+(b23_02+b23_01+XI(1)*(b23_xi2+b23_xi1))/2;
b13=(b13_02-b13_01+XI(2)*(b13_eta2-b13_eta1))*XI(3)/2+(b13_02+b13_01+XI(2)*(b13_eta2+b13_eta1))/2;
%% Bc0
     COL=3*(I-1)+1:3*(I-1)+3;
    Bc0(:,COL)=[N1(I)*J1';
                N2(I)*J2';
                (B33A'+B33B'+B33C')/3;
                N1(I)*J2'+N2(I)*J1';
                b23;
                b13];

%     Bc0(:,COL)=[N1(I)*J1';
%                 N2(I)*J2';
%                 N3(I)*J3';
%                 N1(I)*J2'+N2(I)*J1';
%                 N3(I)*J2'+N2(I)*J3';
%                 N1(I)*J3'+N3(I)*J1'];
end
end