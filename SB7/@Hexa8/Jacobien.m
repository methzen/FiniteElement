function Jac=Jacobien(obj,XYZ,Xsi,Eta,Zeta)
% d�rive X, Y ,Z par rapport Xsi Eta Zeta
X=XYZ(:,1);
Y=XYZ(:,2);
Z=XYZ(:,3);
[dN_dXsi,dN_dEta,dN_dZeta]=DShape8(obj,Xsi,Eta,Zeta);

dX_dXsi =dot(X,dN_dXsi);
dX_dEta =dot(X,dN_dEta);
dX_dZeta=dot(X,dN_dZeta);

dY_dXsi =dot(Y,dN_dXsi);
dY_dEta =dot(Y,dN_dEta);
dY_dZeta=dot(Y,dN_dZeta);

dZ_dXsi =dot(Z,dN_dXsi);
dZ_dEta =dot(Z,dN_dEta);
dZ_dZeta=dot(Z,dN_dZeta);


Jac=[
    dX_dXsi     dY_dXsi     dZ_dXsi
    dX_dEta     dY_dEta     dZ_dEta
    dX_dZeta    dY_dZeta    dZ_dZeta
];


end