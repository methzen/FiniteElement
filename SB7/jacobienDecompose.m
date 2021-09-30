function [J10,J1ZETA,...
          J20,J2ZETA,...
          J30,J3ETA,J3XI]=jacobienDecompose(XYZ)

% XYZ == Coordonnées initiales
%
U=XYZ(:,1);
V=XYZ(:,2);
W=XYZ(:,3);
%
[g1,g2,g3,h1,h2]=g1g2g3h1h2;
%
J10=[dot(g1,U) dot(g1,V) dot(g1,W)];
J1ZETA=[dot(h2,U) dot(h2,V) dot(h2,W)];
%
J20=[dot(g2,U) dot(g2,V) dot(g2,W)];
J2ZETA=[dot(h1,U) dot(h1,V) dot(h1,W)];
%
J30=[dot(g3,U) dot(g3,V) dot(g3,W)];
J3XI=[dot(h2,U) dot(h2,V) dot(h2,W)];
J3ETA=[dot(h1,U) dot(h1,V) dot(h1,W)];
%
end