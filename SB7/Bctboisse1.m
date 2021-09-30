function Bct=Bctboisse1(XYZ,XI)
% face sup
Bctsup=Bctboisse(XYZ,1);
% face inf
Bctinf=Bctboisse(XYZ,-1);
% Interpolation dans l'epaisseur
Bct=(0.5*(1-XI(3))*Bctinf+0.5*(1+XI(3))*Bctsup);
%Bct=(5/6)*(1-XI(3)*XI(3))*(-Bctsup+Bctinf)/2;
end