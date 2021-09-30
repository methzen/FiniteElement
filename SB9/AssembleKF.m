function [GKF,FORCE]=AssembleKF(PROP,coef,LTAN, XYZ,DISPTD,FREEDOF,Connec)
%***********************************************************************
%  MAIN PROGRAM COMPUTING GLOBAL STIFFNESS MATRIX RESIDUAL FORCE FOR
%  Hyperelastic material model assembly
%***********************************************************************
%  initialisation 
   [NUMNP, NDOF] = size(XYZ);				% Analysis parameters
   NE = size(Connec,1);
   p1=NUMNP+1;
   p2=NUMNP+NE;
   DDLsup=p1:p2;
   NEQ = NDOF*NUMNP;
   LE=[Connec DDLsup'];
   GKF = sparse(NEQ+NE,NEQ+NE); % the tangent matrix
   FORCE = sparse(NEQ+NE,1); % RESIDUAL FORCE
   % Geometrie actuelle
   DISPTD=coef*DISPTD;
   xyz=XYZ;
   xyz(FREEDOF) = XYZ(FREEDOF) + DISPTD(FREEDOF);
  %
  % Integration points and weights
  NPG=5;
  [ZG,WGT]=lobatto(NPG);
  %ZG=[-1 -sqrt(3/7) 0 sqrt(3/7)  1];
  %WGT=[0.1 49.0/90.0 32.0/45.0 49.0/90.0 0.1];
  %
  % Index for history variables (each integration pt)
  INTN=0;
  %
  %LOOP OVER ELEMENTS, THIS IS MAIN LOOP TO COMPUTE K AND F
  for IE=1:NE
    % Nodal coordinates and incremental displacements
    ELXY=XYZ(LE(IE,1:8),:);
    ELxy=xyz(LE(IE,1:8),:);
    % Local to global mapping
    IDOF=zeros(1,25);
    for I=1:8
      II=(I-1)*NDOF+1;
      IDOF(II:II+2)=(LE(IE,I)-1)*NDOF+1:(LE(IE,I)-1)*NDOF+3;
    end
    %%%%%%%%%%%%%%%%
    IDOF(25)=NEQ+IE;
    %%%%%%%%%%%%%%%%
    DSP=DISPTD(IDOF(1:24));
    w9=DISPTD(IDOF(25));
    DSP=reshape(DSP,NDOF,8);
    %LOOP OVER INTEGRATION POINTS
    Ueff=0.0;
    [B0,BZETA,BZZ,BXIdev,BETAdev,BETAZETAdev,BXIZETAdev]=matriceBcart(ELXY,ELxy);
    for LZ=1:NPG
      E3=ZG(LZ);
      INTN = INTN + 1;
      %
      % Determinant and shape function derivatives
      J0=Jacobien(ELXY,[0 0 0]);
      DET= det(J0);
      FAC=4*WGT(LZ)*DET;
      %
      % deformation de Gren Lagrange centrée
      E0=deformEcart(DSP,ELXY,[0 0 E3]);
      B9=calculB9(ELXY,[0 0 E3]);
      E0=E0+B9*w9;
      % Computer stress and tangent stiffness
      DTAN= matrice_D(PROP(1), PROP(2));
      % Loi de comportement de Saint Venant Kirchhoff
      STRESS=DTAN*E0;
      % Calcul du module de cisaillement effectif
      ueff=calculUeff(STRESS,E0);
      if ueff< 1e-9
          ueff=DTAN(5,5);
      end
      Ueff=Ueff+ueff*WGT(LZ)/2;
      % calcul des matrice de rigidité
      Bm=B0+E3*BZETA+E3*E3*BZZ;
      B=[Bm B9];
      Km=B'*DTAN*B*FAC;
      Kg=rigigeo(STRESS,ELXY,[0 0 E3]);
      Km(1:24,1:24)=Km(1:24,1:24)+FAC*Kg;
      %% pour force interne
      FORCE(IDOF)=FORCE(IDOF)+FAC*B'*STRESS;
      % Tangent stiffness
      if LTAN
        GKF(IDOF,IDOF)=GKF(IDOF,IDOF)+Km;
      end
    end;
     IDOFS=IDOF(1:24);
     Kmstab=DET*rigistab(BXIdev,BETAdev,BETAZETAdev,BXIZETAdev,Ueff);
     FORCE(IDOFS)=FORCE(IDOFS)+Kmstab*DISPTD(IDOFS);
    if LTAN
      Kgstab=DET*rigigeostab(DSP,ELXY,Ueff);
      GKF(IDOFS,IDOFS)=GKF(IDOFS,IDOFS)+Kgstab+Kmstab;
    end
  end


end
