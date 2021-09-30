function Jac=Jacobien(XYZ,XI)
% derive X, Y ,Z par rapport Xsi Eta Zeta
% XI=[Xsi,Eta,Zeta];


% X=XYZ(:,1);
% Y=XYZ(:,2);
% Z=XYZ(:,3);
% 
% DN=dshapePrisme(XI);
% 
% dX_dXsi =dot(X,DN(1,:));
% dX_dEta =dot(X,DN(2,:));
% dX_dZeta=dot(X,DN(3,:));
% 
% dY_dXsi =dot(Y,DN(1,:));
% dY_dEta =dot(Y,DN(2,:));
% dY_dZeta=dot(Y,DN(3,:));
% 
% dZ_dXsi =dot(Z,DN(1,:));
% dZ_dEta =dot(Z,DN(2,:));
% dZ_dZeta=dot(Z,DN(3,:));
% 
% 
% Jac=[
%     dX_dXsi     dY_dXsi     dZ_dXsi
%     dX_dEta     dY_dEta     dZ_dEta
%     dX_dZeta    dY_dZeta    dZ_dZeta
% ];


% 
[J1,J2,J3]=colomnejacobien(XYZ,XI);
Jac=[J1,J2,J3];
Jac=Jac';
end