function E0=deformEcart(DSP,XYZ,XI)

% XYZ == Coordonnées initiales
% DSP == Déplacement total
[Ec0,Eczeta,Eczz,~,~,~,~]=deformatEc(DSP,XYZ);
[T0,~,~,TZETA]=matriceTdecomp(XYZ);

%% DEFORMATION DEPENDANT DE ZETA ET CONSTANTE
E0=T0*Ec0+XI(3)*(T0*Eczeta+TZETA*Ec0)+XI(3)*XI(3)*(T0*Eczz+TZETA*Eczeta);

end