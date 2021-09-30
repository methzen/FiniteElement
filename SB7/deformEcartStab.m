function [EH1,EH2,EH12,EH23,EH13]=deformEcartStab(DSP,XYZ)

%[Ec0,Eczeta,Eczz,EcXI,EcETA,EcETAZETA,EcXIZETA]=deformatEc(DSP,XYZ);
[Ec0,Eczeta,Eczz,EcXI,EcETA,EcXIETA,EcETAZETA,EcXIZETA]=deformatEc(DSP,XYZ);
[T0,TXI,TETA,TZETA]=matriceTdecomp(XYZ);

%% DEFORMATIONS DE STABILISATION
EH1=T0*EcXI+TXI*Ec0;
EH2=T0*EcETA+TETA*Ec0;
EH23=T0*EcETAZETA+TETA*Eczeta+TZETA*EcETA;
EH13=T0*EcXIZETA+TXI*Eczeta+TZETA*EcXI;
EH12=T0*EcXIETA+TXI*EcETA+TETA*EcXI;
end