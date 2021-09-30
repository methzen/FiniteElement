function E=deformEcart(DSP,XYZ,XI)

% XYZ == Coordonnées initiales
% DSP == Déplacement total
Ec=deformatEc(DSP,XYZ,XI);
[E13,E23]=Echassumect(DSP,XYZ,XI);
Ec(5)=E23;
Ec(6)=E13;
T=matriceT(XYZ,XI);
E=T*Ec;
end