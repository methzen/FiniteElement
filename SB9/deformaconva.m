function [D10,D1ETA,D1ZETA,D1ETAZETA,...
          D20,D2XI,D2ZETA,D2XIZETA,...
          D30,D3ETA,D3XI,D3XIETA]=deformaconva(DSP)

U=DSP(1,:);
V=DSP(2,:);
W=DSP(3,:);
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
D10=[dot(g1,U) dot(g1,V) dot(g1,W)];
D1ETA=[dot(h1,U) dot(h1,V) dot(h1,W)];
D1ZETA=[dot(h3,U) dot(h3,V) dot(h3,W)];
D1ETAZETA=[dot(h4,U) dot(h4,V) dot(h4,W)];
%
D20=[dot(g2,U) dot(g2,V) dot(g2,W)];
D2XI=[dot(h1,U) dot(h1,V) dot(h1,W)];
D2ZETA=[dot(h2,U) dot(h2,V) dot(h2,W)];
D2XIZETA=[dot(h4,U) dot(h4,V) dot(h4,W)];
%
D30=[dot(g3,U) dot(g3,V) dot(g3,W)];
D3ETA=[dot(h2,U) dot(h2,V) dot(h2,W)];
D3XI=[dot(h3,U) dot(h3,V) dot(h3,W)];
D3XIETA=[dot(h4,U) dot(h4,V) dot(h4,W)];

end

