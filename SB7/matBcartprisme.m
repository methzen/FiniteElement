function Bpris=matBcartprisme(XYZ,xyz,XI)

% ICI IL FAUT RENTRER AVEC LES POSITIONS ACTUELLES xyz et Initiale XYZ
%
% [Bc0,BcZETA,BcZZ,BcXI,BcETA,BcETAZETA,BcXIZETA]=matriceBc(xyz);
%  Bcpri=Bc0+BcZETA*XI(3)+BcZZ*XI(3)*XI(3)+BcXI*XI(1)+BcETA*XI(2)+...
%        BcETAZETA*XI(2)*XI(3)+BcXIZETA*XI(1)*XI(3);
   
   
Bcpri=matBc(xyz,XI);
T=matriceT(XYZ,XI);
Bpris=T*Bcpri;

end
