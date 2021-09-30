function [SH1,SH2,SH12,SH23,SH13]=contraintestab(DSP,XYZ,C)
% CALCUL DES CONTRAINTES DE STABILISATION 
%
%
%[EH1,EH2,EH23,EH13]=deformEcartStab(DSP,XYZ);
[EH1,EH2,EH12,EH23,EH13]=deformEcartStab(DSP,XYZ);
% On considère la partie déviatorique de ces deformations
%
I=[1 1 1 0 0 0]';
EH1=EH1-(EH1(1)+EH1(2)+EH1(3))*I/3;
EH2=EH2-(EH2(1)+EH2(2)+EH2(3))*I/3;
EH13=EH13-(EH13(1)+EH13(2)+EH13(3))*I/3;
EH23=EH23-(EH23(1)+EH23(2)+EH23(3))*I/3;
EH12=EH12-(EH12(1)+EH12(2)+EH12(3))*I/3;
%

SH1=C*EH1;
SH2=C*EH2;
SH13=C*EH13;
SH23=C*EH23;
SH12=C*EH12;
end