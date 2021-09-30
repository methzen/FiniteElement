function [B0,BZETA,BZZ,BXIdev,BETAdev,BXIETAdev,BETAZETAdev,BXIZETAdev]=matriceBcart(XYZ,xyz)

% ICI IL FAUT RENTRER AVEC LES POSITIONS ACTUELLES xyz et Initiale XYZ
%
[Bc0,BcZETA,BcZZ,BcXI,BcETA,BcXIETA,BcETAZETA,BcXIZETA]=matriceBc(xyz);
[T0,TXI,TETA,TZETA]=matriceTdecomp(XYZ);
%
B0=T0*Bc0;
BZETA=T0*BcZETA+TZETA*Bc0;
BZZ=T0*BcZZ+TZETA*BcZETA;
%
BXI=T0*BcXI+TXI*Bc0;
BETA=T0*BcETA+TETA*Bc0;
BETAZETA=T0*BcETAZETA+TETA*BcZETA+TZETA*BcETA;
BXIETA=T0*BcXIETA+TXI*BcETA+TETA*BcXI;
BXIZETA=T0*BcXIZETA+TXI*BcZETA+TZETA*BcXI;
%
% %% Enlever les parties deviatoriques sur les B de stabilisation
BXIv=zeros(6,24);
%%
for i=1:8
alpha1=0;
alpha2=0;
alpha3=0;    
    col=3*(i-1)+1:3*(i-1)+3;
    B=BXI(1:3,col);
    for j=1:3
        alpha1=alpha1+B(j,1);
        alpha2=alpha2+B(j,2);
        alpha3=alpha3+B(j,3);
    end
    BXIv(:,col)=[alpha1 alpha2 alpha3;
                 alpha1 alpha2 alpha3;
                 alpha1 alpha2 alpha3;
                  0       0     0    ;
                  0       0     0    ;
                  0       0     0     ]/3;
end
BXIdev=BXI-BXIv;
%% 
BETAv=zeros(6,24);
for i=1:8
alpha1=0;
alpha2=0;
alpha3=0;    
    col=3*(i-1)+1:3*(i-1)+3;
    B=BETA(1:3,col);
    for j=1:3
        alpha1=alpha1+B(j,1);
        alpha2=alpha2+B(j,2);
        alpha3=alpha3+B(j,3);
    end
    BETAv(:,col)=[alpha1 alpha2 alpha3;
                  alpha1 alpha2 alpha3;
                  alpha1 alpha2 alpha3;
                  0       0     0    ;
                  0       0     0    ;
                  0       0     0     ]/3;
end
BETAdev=BETA-BETAv;
%%
BXIETAv=zeros(6,24);
for i=1:8
alpha1=0;
alpha2=0;
alpha3=0;   
    col=3*(i-1)+1:3*(i-1)+3;
        B=BXIETA(1:3,col);
    for j=1:3
        alpha1=alpha1+B(j,1);
        alpha2=alpha2+B(j,2);
        alpha3=alpha3+B(j,3);
    end
    BXIETAv(:,col)=[alpha1 alpha2 alpha3;
                      alpha1 alpha2 alpha3;
                      alpha1 alpha2 alpha3;
                       0       0     0    ;
                       0       0     0    ;
                       0       0     0     ]/3;
end
BXIETAdev=BXIETA-BXIETAv;
%%
BETAZETAv=zeros(6,24);
for i=1:8
alpha1=0;
alpha2=0;
alpha3=0;    
    col=3*(i-1)+1:3*(i-1)+3;
        B=BETAZETA(1:3,col);
    for j=1:3
        alpha1=alpha1+B(j,1);
        alpha2=alpha2+B(j,2);
        alpha3=alpha3+B(j,3);
    end
    BETAZETAv(:,col)=[alpha1 alpha2 alpha3;
                      alpha1 alpha2 alpha3;
                      alpha1 alpha2 alpha3;
                       0       0     0    ;
                       0       0     0    ;
                       0       0     0     ]/3;
end
BETAZETAdev=BETAZETA-BETAZETAv;
%
BXIZETAv=zeros(6,24);
  
for i=1:8
alpha1=0;
alpha2=0;
alpha3=0;    
    col=3*(i-1)+1:3*(i-1)+3;
    B=BXIZETA(1:3,col);
    for j=1:3
        alpha1=alpha1+B(j,1);
        alpha2=alpha2+B(j,2);
        alpha3=alpha3+B(j,3);
    end
    BXIZETAv(:,col)=[alpha1 alpha2 alpha3;
                     alpha1 alpha2 alpha3;
                     alpha1 alpha2 alpha3;
                      0       0     0    ;
                      0       0     0    ;
                      0       0     0     ]/3;
end
BXIZETAdev=BXIZETA-BXIZETAv;
end
