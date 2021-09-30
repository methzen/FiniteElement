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
   %LE=Connec;
  % Integration points and weights
   ZG=[-1 -sqrt(3/7) 0 sqrt(3/7)  1];
  WGT=[0.1 49.0/90.0 32.0/45.0 49.0/90.0 0.1];

  %%
  nddl=19;
  % Index for history variables (each integration pt)
  INTN=0;
  %
  %LOOP OVER ELEMENTS, THIS IS MAIN LOOP TO COMPUTE K AND F
  for IE=1:NE
      G7=zeros(3);
    % Nodal coordinates and incremental displacements
    ELXY=XYZ(LE(IE,1:6),:);
    % Local to global mapping
    IDOF=zeros(1,nddl);
    for I=1:6
      II=(I-1)*NDOF+1;
      IDOF(II:II+2)=(LE(IE,I)-1)*NDOF+1:(LE(IE,I)-1)*NDOF+3;
    end
    %%%%%%%%%%%%%%%%
    %IDOF(19:21)=NEQ+3*(IE-1)+1:NEQ+3*(IE-1)+3;
    IDOF(19)=NEQ+IE;
    %%%%%%%%%%%%%%%%
    DSP=DISPTD(IDOF(1:nddl));
    %W=DISPTD(IDOF(19:21));
    ELxy=ELXY;
    %LOOP OVER INTEGRATION POINTS   
    %
    for LZ=1:size(WGT,2);
      %E3=ZG(LZ);
      INTN = INTN + 1;      
      J0=Jacobien(ELXY,[1/3 1/3 ZG(LZ)]);
      DET= det(J0);
      FAC=WGT(LZ)*DET/2;
      BN=calculBcartassume(ELXY,ELxy,[1/3 1/3 ZG(LZ)]);
      %BN=calculmatBSB6(ELXY,[1/3 1/3 ZG(LZ)]);
      %max(max(BN-B))
      %E0=B*DSP+Benh*W;
      E0=BN*DSP;
      % Computer stress and tangent stiffness
      [DTAN,~]= matrice_D(PROP(1), PROP(2));
      % Loi de comportement de Saint Venant Kirchhoff
      STRESS=DTAN*E0;
      % Update plastic variables
      if UPDATE
        SIGMA(:,INTN)=STRESS;
        continue;
      end
      %% Pour force interne
      FORCE(IDOF)=FORCE(IDOF)-FAC*BN'*STRESS;
      % Tangent stiffness
      if LTAN
         GKF(IDOF,IDOF)=GKF(IDOF,IDOF)+BN'*DTAN*BN*FAC;
      end
    end
  end
end

