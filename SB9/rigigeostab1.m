function Kgstab=rigigeostab1(XYZ,SH1,SH2,SH23,SH13,SH12)

%Coordonnées initiale en entrée
[T0,TXI,TETA,TZETA]=matriceTdecomp(XYZ);

%function [Bc0,BcZETA,BcZZ,BcXI,BcETA,BcETAZETA,BcXIZETA]=matriceBc(XYZ)      
%%
g1=[-1 1 1 -1 -1 1 1 -1]/8;
g2=[-1 -1 1 1 -1 -1 1 1]/8;
g3=[-1 -1 -1 -1 1 1 1 1]/8;
%
h1=[1 -1 1 -1 1 -1 1 -1]/8;
h2=[1 1 -1 -1 -1 -1 1 1]/8;
h3=[1 -1 -1 1 -1 1 1 -1]/8;
h4=[-1 1 -1 1 1 -1 1 -1]/8;
IDEN=[1 1 1 0 0 0]';
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
Kgstab=zeros(24,24);
for I=1:8
    for J=1:8
    LIN=3*(I-1)+1:3*(I-1)+3;
    COL=3*(J-1)+1:3*(J-1)+3;
    %% B0
    aux33=0;
    aux23=0;
    aux13=0;
        for AD=1:4
        DSF = DSHAPE(XAD(AD,:));
        aux33=aux33+DSF(3,I)*DSF(3,J);
        end
        for EH=1:4
        DSF = DSHAPE(XEH(EH,:));
        aux23=aux23+(DSF(2,I)*DSF(3,J)+DSF(3,I)*DSF(2,J))/4;
        end
        for JM=1:4
        DSF = DSHAPE(XJM(JM,:));
        aux13=aux13+(DSF(1,I)*DSF(3,J)+DSF(3,I)*DSF(1,J))/4;
        end

          Gc0=[g1(I)*g1(J);
                g2(I)*g2(J);
                aux33;
                g1(I)*g2(J)+g2(I)*g1(J);
                aux23;
                aux13];
    %% BZETA
    aux23=0;
    aux13=0;
    for EH=1:4
       DSF = DSHAPE(XEH(EH,:));
       aux23=aux23+XEH(EH,3)*(DSF(2,I)*DSF(3,J)+DSF(3,I)*DSF(2,J))/4;
    end
    for JM=1:4
       DSF = DSHAPE(XJM(JM,:));
       aux13=aux13+XJM(JM,3)*(DSF(1,I)*DSF(3,J)+DSF(3,I)*DSF(1,J))/4;
    end
    AUX12=h2(I)*g1(J)+h3(I)*g2(J)+g2(I)*h3(J)+g1(I)*h2(J)  ;
           GcZETA=[ h3(I)*g1(J)+g1(I)*h3(J);
                    h2(I)*g2(J)+g2(I)*h2(J);
                    0 ;
                    AUX12;
                    aux23;
                    aux13];
%% BcXI
    AUX33=0;
    aux23=0;
    for AD=1:4
       DSF = DSHAPE(XAD(AD,:));
       AUX33=AUX33+XAD(AD,1)*(DSF(3,I)*DSF(3,J))/4;
    end
    for EH=1:4
       DSF = DSHAPE(XEH(EH,:));
       aux23=aux23+XEH(EH,1)*(DSF(2,I)*DSF(3,J)+DSF(3,I)*DSF(2,J))/4;
    end
                GcXI=[ 0 ;
                    h1(I)*g2(J)+g2(I)*h1(1);
                    AUX33;
                    h1(I)*g1(J)+g1(I)*h1(J);
                    aux23;
                    0];

%% ETA
    AUX33=0;
    aux13=0;
    for AD=1:4
       DSF = DSHAPE(XAD(AD,:));
       AUX33=AUX33+XAD(AD,2)*(DSF(3,I)*DSF(3,J))/4;
    end
    for JM=1:4
       DSF = DSHAPE(XJM(JM,:));
       aux13=aux13+XJM(JM,2)*(DSF(1,I)*DSF(3,J)+DSF(3,I)*DSF(1,J))/4;
    end
               GcETA=[  h1(I)*g1(J)+g1(I)*h1(J);
                     0 ;
                     AUX33;
                     h1(I)*g2(J)+g2(I)*h1(J);                    
                     0 ;
                     aux13];
%% XIETA
    for AD=1:4
       DSF = DSHAPE(XAD(AD,:));
       AUX33=AUX33+XAD(AD,2)*XAD(AD,1)*(DSF(3,I)*DSF(3,J))/4;
    end
           GcXIETA=[  0;
                      0 ;
                      AUX33 ;
                      2*h1(I)*h1(J);                    
                      0    ;
                      0];
%% ETAZETA
    AUX11=h4(I)*g1(J)+h3(I)*h1(J)+h1(I)*h3(J)+g1(I)*h4(J);
    AUX12=h4(I)*g2(J)+h2(I)*h1(J)+h1(I)*h2(J)+g2(I)*h4(J);
    AUX13=0;
    for JM=1:4
       DSF = DSHAPE(XJM(JM,:));
       AUX13=AUX13+XJM(JM,2)*XJM(JM,3)*(DSF(1,I)*DSF(3,J)+DSF(3,I)*DSF(1,J))/4;
    end
           GcETAZETA=[ AUX11;
                    0    ;
                    0   ;
                    AUX12;                    
                    0    ;
                    AUX13];

%% XIZETA
    AUX22=h4(I)*g2(J)+h2(I)*h1(J)+h1(I)*h2(J)+g2(I)*h4(J);
    AUX12=h4(I)*g1(J)+h3(I)*h1(J)+h1(I)*h3(J)+g1(I)*h4(J);
    AUX23=0;
    for EH=1:4
       DSF = DSHAPE(XEH(EH,:));
       AUX23=AUX23+XEH(EH,1)*XEH(EH,3)*(DSF(2,I)*DSF(3,J)+DSF(3,I)*DSF(2,J))/4;
    end
          GcXIZETA=[ 0 ;
                    AUX22;
                    0   ;                    
                    AUX12;
                    AUX23;
                    0   ];
%%   
GXI=T0*GcXI+TXI*Gc0;
GETA=T0*GcETA+TETA*Gc0;
GETAZETA=T0*GcETAZETA+TETA*GcZETA+TZETA*GcETA;
GXIZETA=T0*GcXIZETA+TXI*GcZETA+TZETA*GcXI;
GXIETA=T0*GcXIETA+TXI*GcETA+TETA*GcXI;
% CALCUL DES DEVIATORIQUES
GXI=GXI-(GXI(1)+GXI(2)+GXI(3))*IDEN/3; 
GETA=GETA-(GETA(1)+GETA(2)+GETA(3))*IDEN/3; 
GETAZETA=GETAZETA-(GETAZETA(1)+GETAZETA(2)+GETAZETA(3))*IDEN/3; 
GXIZETA=GXIZETA-(GXIZETA(1)+GXIZETA(2)+GXIZETA(3))*IDEN/3;
GXIETA=GXIETA-(GXIETA(1)+GXIETA(2)+GXIETA(3))*IDEN/3;
% %
CONST=SH1'*GXI*(8/3)+SH2'*GETA*(8/3)+(SH23'*GETAZETA+SH13'*GXIZETA+SH12'*GXIETA)*(8/9);
Kgstab(LIN,COL)=CONST*eye(3);
    end
end
end
