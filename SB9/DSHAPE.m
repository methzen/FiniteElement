function DSF = DSHAPE(XI)
%*************************************************************************
% Compute shape function, derivatives, and determinant of hexahedron element
% XI=[xi,eta,zeta] location in ref coordinates where shape fun are calculated
% ELXY [3x8] contains the nodal coordinates
% SF [1x8] Shape Function
% GDSF [3x8] Derivative shape function
%*************************************************************************
%%
%   XNODE=[-1  1  1 -1 -1  1  1 -1;
%          -1 -1  1  1 -1 -1  1  1;
%          -1 -1 -1 -1  1  1  1  1]; 
  XNODE=[-1  1  1 -1 -1  1 1 -1;
 -1 -1  1  1 -1 -1 1  1;
 -1 -1 -1 -1  1  1 1  1];
  QUAR = 0.125;
  SF=zeros(8,1);
  DSF=zeros(3,8);
  for I=1:8
    XP = XNODE(1,I);
    YP = XNODE(2,I);
    ZP = XNODE(3,I);
    %
    XI0 = [1+XI(1)*XP 1+XI(2)*YP 1+XI(3)*ZP];
    %
    SF(I) = QUAR*XI0(1)*XI0(2)*XI0(3);
    DSF(1,I) = QUAR*XP*XI0(2)*XI0(3);
    DSF(2,I) = QUAR*YP*XI0(1)*XI0(3);
    DSF(3,I) = QUAR*ZP*XI0(1)*XI0(2);
  end
  
end