function Kg=rigigeo(S,XYZ,XI)

zeta=XI(3);
%Coordonnées initiale en entrée
% Cette matrice Kg est celle sur un points de Gauss
[T0,TXI,TETA,TZETA]=matriceTdecomp(XYZ);
%function [Bc0,BcZETA,BcZZ,BcXI,BcETA,BcETAZETA,BcXIZETA]=matriceBc(XYZ)
%%
g1=[-1  1  1 -1 -1  1 1 -1]/8;
g2=[-1 -1  1  1 -1 -1 1  1]/8;
g3=[-1 -1 -1 -1  1  1 1  1]/8;
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
Kg=zeros(24,24);
for I=1:8
    for J=1:8
    LIN=3*(I-1)+1:3*(I-1)+3;
    COL=3*(J-1)+1:3*(J-1)+3;

    aux033=0;
    aux023=0;
    aux013=0;
    auxZ23=0;
    auxZ13=0;
   %% B0
        for AD=1:4
        DSF = DSHAPE(XAD(AD,:));
        aux033=aux033+DSF(3,I)*DSF(3,J)/4;
%
        DSF = DSHAPE(XEH(AD,:));
        aux023=aux023+(DSF(2,I)*DSF(3,J)+DSF(3,I)*DSF(2,J))/4;
        auxZ23=auxZ23+XEH(AD,3)*(DSF(2,I)*DSF(3,J)+DSF(3,I)*DSF(2,J))/4;
%
        DSF = DSHAPE(XJM(AD,:));
        aux013=aux013+(DSF(1,I)*DSF(3,J)+DSF(3,I)*DSF(1,J))/4;
        auxZ13=auxZ13+XJM(AD,3)*(DSF(1,I)*DSF(3,J)+DSF(3,I)*DSF(1,J))/4;
        end
        %mm
        %aux033
           GIJ0=[g1(I)*g1(J);
                 g2(I)*g2(J);
                 aux033;
                 g1(I)*g2(J)+g2(I)*g1(J);
                 aux023;
                 aux013];        
    %% BZETA
       AUXZ12=h2(I)*g1(J)+h3(I)*g2(J)+g2(I)*h3(J)+g1(I)*h2(J)  ;
              GIJZ=[h3(I)*g1(J)+g1(I)*h3(J);
                    h2(I)*g2(J)+g2(I)*h2(J);
                    0 ;
                    AUXZ12;
                    auxZ23;
                    auxZ13];
    %% BZZ
             GIJZZ=[h3(I)*h3(J);
                    h2(I)*h2(J);
                    0 ;
                    h2(I)*h3(J)+h3(I)*h2(J);
                    0 ;
                    0 ];
%%
        alpha=T0*GIJ0+zeta*(T0*GIJZ+TZETA*GIJ0)+zeta*zeta*(T0*GIJZZ+TZETA*GIJZ); 
        Kg(LIN,COL)=dot(alpha,S)*eye(3);
    end
end

end
