function Jac=Jacobien(XYZ,XI)
% derive X, Y ,Z par rapport Xsi Eta Zeta
% XI=[Xsi,Eta,Zeta];
X=XYZ(:,1);
Y=XYZ(:,2);
Z=XYZ(:,3);
g1=[-1 1 1 -1 -1 1 1 -1]/8;
g2=[-1 -1 1 1 -1 -1 1 1]/8;
g3=[-1 -1 -1 -1 1 1 1 1]/8;
%
h1=[1 -1 1 -1 1 -1 1 -1]/8;
h2=[1 1 -1 -1 -1 -1 1 1]/8;
h3=[1 -1 -1 1 -1 1 1 -1]/8;
h4=[-1 1 -1 1 1 -1 1 -1]/8;

%
dN_dXsi =g1+XI(2)*h1+XI(3)*h3+XI(2)*XI(3)*h4;
dN_dEta =g2+XI(1)*h1+XI(3)*h2+XI(1)*XI(3)*h4;
dN_dZeta=g3+XI(2)*h2+XI(1)*h3+XI(2)*XI(1)*h4;
%
dX_dXsi =dot(X,dN_dXsi);
dX_dEta =dot(X,dN_dEta);
dX_dZeta=dot(X,dN_dZeta);
%
dY_dXsi =dot(Y,dN_dXsi);
dY_dEta =dot(Y,dN_dEta);
dY_dZeta=dot(Y,dN_dZeta);
%
dZ_dXsi =dot(Z,dN_dXsi);
dZ_dEta =dot(Z,dN_dEta);
dZ_dZeta=dot(Z,dN_dZeta);
%
Jac=[
    dX_dXsi     dY_dXsi     dZ_dXsi
    dX_dEta     dY_dEta     dZ_dEta
    dX_dZeta    dY_dZeta    dZ_dZeta
];
Jac=Jac';
end