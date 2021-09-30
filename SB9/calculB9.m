function B9=calculB9(XYZ,XI)
% Calcul de la composante du neuvieme noeuds
[T0,~,~,~]=matriceTdecomp(XYZ);
Bc9=[0 0 XI(3) 0 0 0]';
B9=T0*Bc9;
end