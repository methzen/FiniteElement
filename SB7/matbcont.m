function B=matbcont(XYZ,XI)
[Bc0,BcZETA,BcZZ,BcXI,BcETA,BcETAZETA,BcXIZETA]=matriceBc(XYZ);

B=Bc0+BcZETA*XI(3)+BcZZ*XI(3)*XI(3)+BcXI*XI(1)+BcETA*XI(2)+BcETAZETA*XI(2)*XI(3)+BcXIZETA*XI(1)*XI(3);
end