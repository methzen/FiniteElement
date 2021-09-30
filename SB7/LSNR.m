function LSNR(CON,FORCETOT,SDISPT,PROP,XYZ,LE)

global DISPTD
[NUMNP, NDOF] = size(XYZ);
 NEQ = NDOF*NUMNP+size(LE,1);
FIXEDDOF=3*(SDISPT(:,1)-1)+SDISPT(:,2);
ALLDOF=1:NEQ;
FREEDOF=setdiff(ALLDOF,FIXEDDOF);
FORCEXTG
% Starts the load increment loop.
%--------------------------------------------------------------------------
while ((CON.xlamb < CON.xlmax) && (CON.incrm < CON.nincr))
      %--------------------------------------------------------------------
      % Update the increment.
      %--------------------------------------------------------------------
      CON.incrm = CON.incrm + 1;
      %--------------------------------------------------------------------
      % Update the load factor. 
      %--------------------------------------------------------------------
      CON.xlamb = CON.xlamb + CON.dlamb;
      %--------------------------------------------------------------------
      % Update nodal forces (excluding pressure) and gravity. 
      %--------------------------------------------------------------------
      [RESIDU,FORCEXTG] = external_force_update(FORCETOT,RESIDU,FORCEXTG,CON.dlamb);
      %--------------------------------------------------------------------
      [GKF,~]=AssembleKF(PROP,1,LTAN, XYZ,DISPTD,FREEDOF,LE);
      % Newton-Raphson iteration.     
      %--------------------------------------------------------------------
      CON.niter = 0;
      rnorm = 2*CON.cnorm;
      while((rnorm > CON.cnorm) && (CON.niter < CON.miter))
          CON.niter = CON.niter + 1;
          %----------------------------------------------------------------
          % Solve for iterative displacements. Also obtains the product r.u.
          %----------------------------------------------------------------
          [SOL,rtu0] = linear_solver(GKF,-RESIDU,BC.fixdof);
          %----------------------------------------------------------------
          % Starts the line search iteration. The total number of line search
          % iterations is limited to msearch.  
          %----------------------------------------------------------------
          eta0  = 0;
          eta   = 1;   
          nsear = 0;  
          rtu   = rtu0*CON.searc*2;
          while((abs(rtu)  > abs(rtu0*CON.searc) ) && (nsear < CON.msearch))
              nsear = nsear + 1; 
              %------------------------------------------------------------
              % Update coodinates.
              % -Recompute equivalent nodal forces and assembles residual 
              %  force, excluding pressure contributions.
              % -Recompute and assembles tangent stiffness matrix components,
              %  excluding pressure contributions.
              %------------------------------------------------------------
              alpha=eta-eta0;
              [~,FORCE]=AssembleKF(PROP,alpha,LTAN, XYZ,DISPTD,FREEDOF,Connec);
              % calcul de Fl=xlamb*F
              % calcul de GKF avec geom a jour
              % Calcul de Rl=Fint-FL
              RESIDU=xlamb*FORCEXTG-FORCE;
              %------------------------------------------------------------
              % Line-search algorithm.
              %------------------------------------------------------------
              [eta0,eta,rtu] = search(eta0,eta,rtu0,RESIDU, SOL, FREEDOF);
          end     
          %----------------------------------------------------------------
          % Check for equilibrium convergence.
          %----------------------------------------------------------------
          [rnorm,GLOBAL] = check_residual_norm(CON,BC,GLOBAL,BC.freedof);
          %----------------------------------------------------------------
          % Break iteration before residual gets unrealistic (e.g. NaN).
          %----------------------------------------------------------------
          if (abs(rnorm)>1e7 || isnan(rnorm))
              CON.niter=CON.miter;
              break;
          end         
      end      
      %--------------------------------------------------------------------
      % If convergence not achieved opt to restart from previous results or
      % terminate.
      %--------------------------------------------------------------------
      if( CON.niter >= CON.miter)          
        terminate = input(['Solution not converged. Do you want to '...
                           'terminate the program (y/n) ?: ' ' \n'],'s');
        if( strcmp(deblank(terminate),deblank('n')) || strcmp(deblank(terminate),deblank('N')))
            load(PRO.internal_restartfile_name);   
            fprintf(['Restart from previous step by decreasing '...
                     'the load parameter increment to half its '...
                     'initally fixed value. \n']);
            CON.dlamb = CON.dlamb/2;
        else
            error('Program terminated by user');
        end
      else
        %------------------------------------------------------------------         
      end
end       
fprintf(' Normal end of PROGRAM flagshyp. \n');


%% affichage
NDOF=3*size(XYZ,1);
DSP=full(DISPTD(1:NDOF));
DSP=reshape(DSP,3,size(XYZ,1));
XYZ2=XYZ+10*DSP';
for i=1:size(LE,1)
plotHexa8(XYZ(LE(i,:),:),'b')
plotHexa8(XYZ2(LE(i,:),:),'r')
end

end