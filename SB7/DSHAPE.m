function DSF = DSHAPE(XI)
%*************************************************************************
% Compute shape function, derivatives, and determinant of hexahedron element
% XI=[xi,eta,zeta] location in ref coordinates where shape fun are calculated
% ELXY [3x8] contains the nodal coordinates
% SF [1x8] Shape Function
% GDSF [3x8] Derivative shape function
%*************************************************************************
%%
[g1,g2,g3,h1,h2]=g1g2g3h1h2;

DSF(1,:) = g1+h1*XI(3);
DSF(2,:) = g2+h2*XI(3);
DSF(3,:) = g3+h1*XI(1)+h2*XI(2);
  
end