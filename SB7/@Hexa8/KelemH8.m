function Ke=KelemH8(obj,X,Y,Z,E,v)
% Calcul de la rigidité élémentaire en fonction de Lamda et mu

% initialisation
Ke=sparse(24,24); 

% position des points de Gauss
XsiG  =1./sqrt(3)*[-1  1  1 -1 -1  1  1 -1];
EtaG  =1./sqrt(3)*[-1 -1  1  1 -1 -1  1  1];
ZetaG =1./sqrt(3)*[-1 -1 -1 -1  1  1  1  1];
% poids de Gauss
weightGauss=1;

%% boucle sur les points de Gauss pour intégration complète
for i=1:8
    % calcul du determinent du Jacobien
    DetJac=det(Jacobien(obj,X,Y,Z,XsiG(i),EtaG(i),ZetaG(i)));
    
    % Calcul de la matrice de rigidité élémentaire
    Ke=Ke+...
    weightGauss*DetJac*BepsTxDxBeps(obj,X,Y,Z,XsiG(i),EtaG(i),ZetaG(i),E,v);

end

%%
function BtDB=BepsTxDxBeps(obj,X,Y,Z,Xsi,Eta,Zeta,E,v)
% calcul de la Bt * D* B en Xsi, Eta, Zeta

% matrice de comportement
lambda=v*E/((1-2*v)*(1+v));
mu=E/(2*(1+v));
D=[
    (lambda+2*mu)   lambda          lambda          0       0       0
    lambda          (lambda+2*mu)   lambda          0       0       0
    lambda          lambda          (lambda+2*mu)   0       0       0
    0               0               0               mu      0       0
    0               0               0               0       mu      0
    0               0               0               0       0       mu  
];

% matrice B
Beps=calculeB(obj,X,Y,Z,Xsi,Eta,Zeta);

% matrice Kg
BtDB=Beps'*D*Beps;

end

end
