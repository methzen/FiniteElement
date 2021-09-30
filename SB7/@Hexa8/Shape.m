function N=Shape(obj,Xsi, Eta, Zeta)
% Calcul des 8 fonctions de formes 

N(1)=0.125*(1-Xsi)*(1-Eta)*(1-Zeta);
N(2)=0.125*(1+Xsi)*(1-Eta)*(1-Zeta);
N(3)=0.125*(1+Xsi)*(1+Eta)*(1-Zeta);
N(4)=0.125*(1-Xsi)*(1+Eta)*(1-Zeta);
N(5)=0.125*(1-Xsi)*(1-Eta)*(1+Zeta);
N(6)=0.125*(1+Xsi)*(1-Eta)*(1+Zeta);
N(7)=0.125*(1+Xsi)*(1+Eta)*(1+Zeta);
N(8)=0.125*(1-Xsi)*(1+Eta)*(1+Zeta);

end