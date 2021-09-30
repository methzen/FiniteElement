function ELAST3D(PROP,UPDATE, LTAN, XYZ,Connec)
%***********************************************************************
%  MAIN PROGRAM COMPUTING GLOBAL STIFFNESS MATRIX RESIDUAL FORCE FOR
%  Hyperelastic material model assembly
%***********************************************************************
%  initialisation 
global DISPTD FORCE GKF SIGMA
   [NUMNP, NDOF] = size(XYZ);				% Analysis parameters
   NE = size(Connec,1);
   p1=NUMNP+1;
   p2=NUMNP+NE;
   DDLsup=p1:p2;
   NEQ = NDOF*NUMNP;
   LE=[Connec DDLsup'];
  %
  % Integration points and weights
  ZG=[-1 -sqrt(3/7) 0 sqrt(3/7)  1];
  WGT=[0.1 49.0/90.0 32.0/45.0 49.0/90.0 0.1];
  %% matrice tangente de stabilisation(Deviatorique)
  Cstab=mattdev;
  %%
  %
  % Index for history variables (each integration pt)
  INTN=0;
  %
  %LOOP OVER ELEMENTS, THIS IS MAIN LOOP TO COMPUTE K AND F
  for IE=1:NE
    % Nodal coordinates and incremental displacements
    ELXY=XYZ(LE(IE,1:8),:);
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
    ELxy=ELXY;
    %LOOP OVER INTEGRATION POINTS
    Ueff=0.0;    
    [B0,BZETA,BZZ,BXIdev,BETAdev,BXIETAdev,BETAZETAdev,BXIZETAdev]=matriceBcart(ELXY,ELxy);
      % Determinant and shape function derivatives
      J0=Jacobien(ELXY,[0 0 0]);
      DET= det(J0);
    for LZ=1:5
      E3=ZG(LZ);
      INTN = INTN + 1;
      %
      FAC=4*WGT(LZ)*DET;
      %
      % deformation de Gren Lagrange centrée
      Bm=B0+E3*BZETA+E3*E3*BZZ;
      E0=Bm*DSP;
      B9=calculB9(ELXY,[0 0 E3]);
      E0=E0+B9*w9;
      % Computer stress and tangent stiffness
      DTAN= matrice_D(PROP(1), PROP(2));
      % Loi de comportement de Saint Venant Kirchhoff
      STRESS=DTAN*E0;
      % Update plastic variables
      if UPDATE
        SIGMA(:,INTN)=STRESS;
        continue;
      end
      % Calcul du module de cisaillement effectif
      Ueff=DTAN(5,5);
      % calcul des matrices de gradient BN et BG
      BN=[Bm B9];     
      %% pour force interne
      FORCE(IDOF)=FORCE(IDOF)-FAC*BN'*STRESS;
      % Tangent stiffness
      if LTAN
         Km=BN'*DTAN*BN*FAC;
                 A=zeros(25,1);
        for i=1:25
            A(i)=Km(i,i);
        end
         GKF(IDOF,IDOF)=GKF(IDOF,IDOF)+Km;
      end
    end;
    if UPDATE==false
     IDOFS=IDOF(1:24);
     Fintstab=DET*(BXIdev'*Cstab*BXIdev/3+BETAdev'*Cstab*BETAdev/3+...
                   BETAZETAdev'*Cstab*BETAZETAdev/9+...
                   BXIZETAdev'*Cstab*BXIZETAdev/9+...
                   0*BXIETAdev'*Cstab*BXIETAdev/9)*Ueff*8;     
     FORCE(IDOFS)=FORCE(IDOFS)-Fintstab*DSP;
    end
    if LTAN
      Kmstab=rigistab(BXIdev,BETAdev,BXIETAdev*0,BETAZETAdev,BXIZETAdev,Ueff*Cstab);
      GKF(IDOFS,IDOFS)=GKF(IDOFS,IDOFS)+DET*Kmstab;
    end
  end

end

