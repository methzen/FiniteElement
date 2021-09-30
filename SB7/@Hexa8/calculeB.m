function Beps=calculeB(obj,XYZ,Xsi,Eta,Zeta)
X=XYZ(:,1);
Y=XYZ(:,2);
Z=XYZ(:,3);

zero8=zeros(1,8);
[dN_dXsi,dN_dEta,dN_dZeta]=DShape8(obj,Xsi,Eta,Zeta);
Jac=Jacobien(obj,XYZ,Xsi,Eta,Zeta);
Jac_Inv=inv(Jac);

%% D�riv�e des coordonn�es isoparam�triques par � X,Y,Z

dXsi_DX=Jac_Inv(1,1);
dXsi_DY=Jac_Inv(2,1);
dXsi_DZ=Jac_Inv(3,1);

dEta_DX=Jac_Inv(1,2);
dEta_DY=Jac_Inv(2,2);
dEta_DZ=Jac_Inv(3,2);

dZeta_DX=Jac_Inv(1,3);
dZeta_DY=Jac_Inv(2,3);
dZeta_DZ=Jac_Inv(3,3);

%% calcul des d�riv�es de N par rapport � X,Y,Z
dN_dX= dN_dXsi*dXsi_DX+dN_dEta*dEta_DX+dN_dZeta*dZeta_DX;
dN_dY= dN_dXsi*dXsi_DY+dN_dEta*dEta_DY+dN_dZeta*dZeta_DY;
dN_dZ= dN_dXsi*dXsi_DZ+dN_dEta*dEta_DZ+dN_dZeta*dZeta_DZ;


%% calcul de Beps
Beps=zeros(6,24); %initialisation

Beps(1,:)=[dN_dX , zero8, zero8]; % pour epsilon XX
Beps(2,:)=[zero8 , dN_dY, zero8]; % pour epsilon YY
Beps(3,:)=[zero8 , zero8, dN_dZ]; % pour epsilon ZZ
Beps(4,:)=[dN_dY , dN_dX, zero8]; % pour 2 epsilon XY
Beps(5,:)=[dN_dZ , zero8, dN_dX]; % pour 2 epsilon XZ
Beps(6,:)=[zero8 , dN_dZ, dN_dY]; % pour 2 epsilon YZ


end