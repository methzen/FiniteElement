function [J10,J1ETA,J1ZETA,J1ETAZETA,...
          J20,J2XI, J2ZETA,J2XIZETA,...
          J30,J3ETA,J3XI,  J3XIETA]=jacobienDecompose(XYZ)

% XYZ == Coordonnées initiales
%
U=XYZ(:,1);
V=XYZ(:,2);
W=XYZ(:,3);
%
g1=[-1 1 1 -1 -1 1 1 -1]/8;
g2=[-1 -1 1 1 -1 -1 1 1]/8;
g3=[-1 -1 -1 -1 1 1 1 1]/8;
%
h1=[1 -1 1 -1 1 -1 1 -1]/8;
h2=[1 1 -1 -1 -1 -1 1 1]/8;
h3=[1 -1 -1 1 -1 1 1 -1]/8;
h4=[-1 1 -1 1 1 -1 1 -1]/8;
%
J10=[dot(g1,U) dot(g1,V) dot(g1,W)];
J1ETA=[dot(h1,U) dot(h1,V) dot(h1,W)];
J1ZETA=[dot(h3,U) dot(h3,V) dot(h3,W)];
J1ETAZETA=[dot(h4,U) dot(h4,V) dot(h4,W)];
%
J20=[dot(g2,U) dot(g2,V) dot(g2,W)];
J2XI=[dot(h1,U) dot(h1,V) dot(h1,W)];
J2ZETA=[dot(h2,U) dot(h2,V) dot(h2,W)];
J2XIZETA=[dot(h4,U) dot(h4,V) dot(h4,W)];
%
J30=[dot(g3,U) dot(g3,V) dot(g3,W)];
J3XI=[dot(h3,U) dot(h3,V) dot(h3,W)];
J3ETA=[dot(h2,U) dot(h2,V) dot(h2,W)];
J3XIETA=[dot(h4,U) dot(h4,V) dot(h4,W)];
end