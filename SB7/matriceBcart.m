function [Bc,Br,Bs,Bz,Bzz,Brz,Bsz,Brs]=matriceBcart(XYZ,xyz)

% ICI IL FAUT RENTRER AVEC LES POSITIONS ACTUELLES xyz et Initiale XYZ
%
 [T0,TXI,TETA,TZETA]=matriceTdecomp(XYZ);
 [Bc0,Bczeta,Bczz,BcXI,BcETA,BczetaXI,BczetaETA,BcXIXI,BcETET,BcXIETA]=matriceBc(xyz);
%
Bc=T0*Bc0;
Br=T0*BcXI+TXI*Bc0;
Bs=T0*BcETA+TETA*Bc0;
Bz=T0*Bczeta+TZETA*Bc0;
Bzz=T0*Bczz+TZETA*Bczeta;
Brz=T0*BczetaXI+TZETA*BcXI+TXI*Bczeta;
Bsz=T0*BczetaETA+TZETA*BcETA+TETA*Bczeta;
Brs=T0*BcXIETA+TETA*BcXI+TXI*BcETA;
%Brr=TXI*BcXI+T0*BcXIXI;
%Bss=TETA*BcETA+T0*BcETET;
%Brzz=TZETA*BczetaXI+TXI*Bczz;
%Bszz=TZETA*BczetaETA+TETA*Bczz;
%Bssz=TZETA*BcETET+TETA*BczetaETA;
%Brrz=TZETA*BcXIXI+TXI*BczetaXI;
%Brsz=TZETA*BcXIETA+TXI*BczetaETA+TETA*BczetaXI;
end
