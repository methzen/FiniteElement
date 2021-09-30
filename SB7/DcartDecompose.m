function [D10,D1ZETA,...
          D20,D2ZETA,...
          D30,D3ETA,D3XI]=DcartDecompose(DEP)

% XYZ == Coordonn�es initiales
%
U=DEP(:,1);
V=DEP(:,2);
W=DEP(:,3);
%
[g1,g2,g3,h1,h2]=g1g2g3h1h2;
%
D10=[dot(g1,U) dot(g1,V) dot(g1,W)];
D1ZETA=[dot(h1,U) dot(h1,V) dot(h1,W)];
%
D20=[dot(g2,U) dot(g2,V) dot(g2,W)];
D2ZETA=[dot(h2,U) dot(h2,V) dot(h2,W)];
%
D30=[dot(g3,U) dot(g3,V) dot(g3,W)];
D3XI=[dot(h1,U) dot(h1,V) dot(h1,W)];
D3ETA=[dot(h2,U) dot(h2,V) dot(h2,W)];
end