function [Bc0,BcZETA,BcZZ,BcXI,BcETA,BcXIETA,BcETAZETA,BcXIZETA]=matriceBc(XYZ)

% ICI IL FAUT RENTRER AVEC LES POSITIONS ACTUELLES 
[J10,J1ETA,J1ZETA,J1ETAZETA,...
 J20,J2XI,J2ZETA,J2XIZETA,...
 ~,~,~,~]=jacobienDecompose(XYZ);
      
%%
g1=[-1 1 1 -1 -1 1 1 -1]/8;
g2=[-1 -1 1 1 -1 -1 1 1]/8;
g3=[-1 -1 -1 -1 1 1 1 1]/8;
%
h1=[1 -1 1 -1 1 -1 1 -1]/8;
h2=[1 1 -1 -1 -1 -1 1 1]/8;
h3=[1 -1 -1 1 -1 1 1 -1]/8;
h4=[-1 1 -1 1 1 -1 1 -1]/8;
%
%% TYING POINTS
% A à D
XAD=[-1 -1 0
      1 -1 0
      1  1 0
     -1  1 0];
% E TO H
XEH=[-1  0 -1 
      1  0 -1
      1  0  1
     -1  0  1];
% J TO M
XJM=[0 -1 -1
     0  1 -1
     0  1  1
     0 -1  1];
%%
 Bc0=zeros(6,24); 
for I=1:8
    aux33=zeros(3,1);
    aux23=zeros(3,1);
    aux13=zeros(3,1);
    for AD=1:4
      [~,~,J3]=colomnejacobien(XYZ,XAD(AD,:));
       DSF = DSHAPE(XAD(AD,:));
       aux33=aux33+DSF(3,I)*J3/4;
    end
    for EH=1:4
      [~,J2,J3]=colomnejacobien(XYZ,XEH(EH,:));
       DSF = DSHAPE(XEH(EH,:));
       aux23=aux23+(DSF(2,I)*J3+DSF(3,I)*J2)/4;
    end
    for JM=1:4
      [J1,~,J3]=colomnejacobien(XYZ,XJM(JM,:));
       DSF = DSHAPE(XJM(JM,:));
       aux13=aux13+(DSF(1,I)*J3+DSF(3,I)*J1)/4;
    end
    COL=3*(I-1)+1:3*(I-1)+3;
    Bc0(:,COL)=[g1(I)*J10(1)  g1(I)*J10(2)  g1(I)*J10(3);
                g2(I)*J20(1)  g2(I)*J20(2)  g2(I)*J20(3);
                    aux33(1)      aux33(2)     aux33(3);
                g1(I)*J20(1)+g2(I)*J10(1) g1(I)*J20(2)+g2(I)*J10(2) g1(I)*J20(3)+g2(I)*J10(3);
                aux23(1)      aux23(2)     aux23(3);
                aux13(1)      aux13(2)     aux13(3)];
end
%% BcZETA
BcZETA=zeros(6,24);
for I=1:8
    aux23=zeros(3,1);
    aux13=zeros(3,1);
    for EH=1:4
      [~,J2,J3]=colomnejacobien(XYZ,XEH(EH,:));
       DSF = DSHAPE(XEH(EH,:));
       aux23=aux23+XEH(EH,3)*(DSF(2,I)*J3+DSF(3,I)*J2)/4;
    end
    for JM=1:4
      [J1,~,J3]=colomnejacobien(XYZ,XJM(JM,:));
       DSF = DSHAPE(XJM(JM,:));
       aux13=aux13+XJM(JM,3)*(DSF(1,I)*J3+DSF(3,I)*J1)/4;
    end
    COL=3*(I-1)+1:3*(I-1)+3;
    AUX12=h2(I)*J10+h3(I)*J20+g2(I)*J1ZETA+g1(I)*J2ZETA  ;
    BcZETA(:,COL)=[ h3(I)*J10(1)+g1(I)*J1ZETA(1)  h3(I)*J10(2)+g1(I)*J1ZETA(2)  h3(I)*J10(3)+g1(I)*J1ZETA(3);
                    h2(I)*J20(1)+g2(I)*J2ZETA(1)  h2(I)*J20(2)+g2(I)*J2ZETA(2)  h2(I)*J20(3)+g2(I)*J2ZETA(3);
                    0      0     0;
                    AUX12(1) AUX12(2) AUX12(3);
                    aux23(1)      aux23(2)     aux23(3);
                    aux13(1)      aux13(2)     aux13(3)];
end
%% ZETA ZETA
BcZZ=zeros(6,24);
for I=1:8
    COL=3*(I-1)+1:3*(I-1)+3;
    BcZZ(:,COL)=[h3(I)*J1ZETA(1)  h3(I)*J1ZETA(2)  h3(I)*J1ZETA(3);
                 h2(I)*J2ZETA(1)  h2(I)*J2ZETA(2)  h2(I)*J2ZETA(3);
                 0      0     0;
                 h2(I)*J1ZETA(1)+h3(I)*J2ZETA(1) h2(I)*J1ZETA(2)+h3(I)*J2ZETA(2) h2(I)*J1ZETA(3)+h3(I)*J2ZETA(3);
                 0      0     0;
                 0      0     0];
end
%% BcXI
BcXI=zeros(6,24);
for I=1:8
    AUX33=zeros(3,1);
    aux23=zeros(3,1);
    for AD=1:4
      [~,~,J3]=colomnejacobien(XYZ,XAD(AD,:));
       DSF = DSHAPE(XAD(AD,:));
       AUX33=AUX33+XAD(AD,1)*(DSF(3,I)*J3)/4;
    end
    for EH=1:4
      [~,J2,J3]=colomnejacobien(XYZ,XEH(EH,:));
       DSF = DSHAPE(XEH(EH,:));
       aux23=aux23+XEH(EH,1)*(DSF(2,I)*J3+DSF(3,I)*J2)/4;
    end
    COL=3*(I-1)+1:3*(I-1)+3;
    BcXI(:,COL)=[ 0      0     0;
                    h1(I)*J20(1)+g2(I)*J2XI(1)  h1(I)*J20(2)+g2(I)*J2XI(2)  h1(I)*J20(3)+g2(I)*J2XI(3);
                    AUX33(1)      AUX33(2)     AUX33(3);
                    h1(I)*J10(1)+g1(I)*J2XI(1)  h1(I)*J10(2)+g1(I)*J2XI(2)  h1(I)*J10(3)+g1(I)*J2XI(3);
                    aux23(1)      aux23(2)     aux23(3);
                    0      0     0];
end
%% ETA
BcETA=zeros(6,24);
for I=1:8
    AUX33=zeros(3,1);
    aux13=zeros(3,1);
    for AD=1:4
      [~,~,J3]=colomnejacobien(XYZ,XAD(AD,:));
       DSF = DSHAPE(XAD(AD,:));
       AUX33=AUX33+XAD(AD,2)*(DSF(3,I)*J3)/4;
    end
    for JM=1:4
      [J1,~,J3]=colomnejacobien(XYZ,XJM(JM,:));
       DSF = DSHAPE(XJM(JM,:));
       aux13=aux13+XJM(JM,2)*(DSF(1,I)*J3+DSF(3,I)*J1)/4;
    end
    COL=3*(I-1)+1:3*(I-1)+3;
    BcETA(:,COL)=[  h1(I)*J10(1)+g1(I)*J1ETA(1)  h1(I)*J10(2)+g1(I)*J1ETA(2)  h1(I)*J10(3)+g1(I)*J1ETA(3);
                    0      0     0;
                    AUX33(1)      AUX33(2)     AUX33(3);
                    h1(I)*J20(1)+g2(I)*J1ETA(1)  h1(I)*J20(2)+g2(I)*J1ETA(2)  h1(I)*J20(3)+g2(I)*J1ETA(3);                    
                    0      0     0;
                    aux13(1)      aux13(2)     aux13(3);];
end
%% XIETA
BcXIETA=zeros(6,24);
for I=1:8
    AUX33=zeros(3,1);
    for AD=1:4
      [~,~,J3]=colomnejacobien(XYZ,XAD(AD,:));
       DSF = DSHAPE(XAD(AD,:));
       AUX33=AUX33+XAD(AD,2)*XAD(AD,1)*(DSF(3,I)*J3)/4;
    end
    COL=3*(I-1)+1:3*(I-1)+3;
    BcXIETA(:,COL)=[0      0     0;
                    0      0     0;
                    AUX33(1)      AUX33(2)     AUX33(3);
                    h1(I)*J2XI(1)+h1(I)*J1ETA(1)  h1(I)*J2XI(2)+h1(I)*J1ETA(2)  h1(I)*J2XI(3)+h1(I)*J1ETA(3);                    
                    0      0     0;
                    0      0     0];
end

%% ETA ZETA
BcETAZETA=zeros(6,24);
for I=1:8
    AUX11=h4(I)*J10+h3(I)*J1ETA+h1(I)*J1ZETA+g1(I)*J1ETAZETA;
    AUX12=h4(I)*J20+h2(I)*J1ETA+h1(I)*J2ZETA+g2(I)*J1ETAZETA;
    AUX13=[0 0 0]';
    for JM=1:4
      [J1,~,J3]=colomnejacobien(XYZ,XJM(JM,:));
       DSF = DSHAPE(XJM(JM,:));
       AUX13=AUX13+XJM(JM,2)*XJM(JM,3)*(DSF(1,I)*J3+DSF(3,I)*J1)/4;
    end
    COL=3*(I-1)+1:3*(I-1)+3;
 BcETAZETA(:,COL)=[ AUX11(1)  AUX11(2)  AUX11(3);
                    0      0     0;
                    0      0     0;
                    AUX12(1)  AUX12(2)  AUX12(3);                    
                    0      0     0;
                    AUX13(1)      AUX13(2)     AUX13(3);];
end
%% XI ZETA
BcXIZETA=zeros(6,24);
for I=1:8
    AUX22=h4(I)*J20+h2(I)*J2XI+h1(I)*J2ZETA+g2(I)*J2XIZETA;
    AUX13=h4(I)*J10+h3(I)*J2XI+h1(I)*J1ZETA+g1(I)*J2XIZETA;
    AUX23=[0 0 0]';
    for EH=1:4
      [~,J2,J3]=colomnejacobien(XYZ,XEH(EH,:));
       DSF = DSHAPE(XEH(EH,:));
       AUX23=AUX23+XEH(EH,1)*XEH(EH,3)*(DSF(2,I)*J3+DSF(3,I)*J2)/4;
    end
    COL=3*(I-1)+1:3*(I-1)+3;
 BcXIZETA(:,COL)=[  0      0     0;
                    AUX22(1)  AUX22(2)  AUX22(3);
                    0      0     0;                    
                    AUX13(1)  AUX13(2)  AUX13(3);
                    AUX23(1)  AUX23(2)  AUX23(3);
                    0      0     0];
end
end