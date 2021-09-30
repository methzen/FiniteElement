function linesearch(PROP,SOLN,FEXT,FREEDOF,XYZ,Connec,rtu0)
global DISPTD DISPDD

% on enregistre les déplacements precedente
    DISPTD0=DISPTD;
    DISPDD0=DISPDD;
    msearch=10;
    searc=0.1;
    DISPTD = DISPTD + 1*SOLN;
    RESIDU=forceinterne(PROP,DISPTD,FEXT,XYZ,Connec);
    rtu=full(dot(SOLN(FREEDOF),RESIDU(FREEDOF)));
    
    if abs(rtu) < abs(rtu0*searc)
        DISPDD = DISPDD0 + SOLN;
        DISPTD = DISPTD0 + SOLN;    
    else
          rho0  = 0;
          rho1   = 1;   
          nsear = 0;  
          while((abs(rtu) >=abs(rtu0*searc) )&& (nsear < msearch))
              
                  alpha=(rho1-rho0)/(rtu-rtu0);
                  rho=rho1-alpha*rtu;
                  DISPTD = DISPTD0 + rho*SOLN;
                  RESIDU=forceinterne(PROP,DISPTD,FEXT, XYZ,Connec);
                  rtu=full(dot(SOLN(FREEDOF),RESIDU(FREEDOF)));
                  
                  rho0=rho1;
                  rho1=rho;
                  nsear = nsear + 1 ;
          end 
        DISPDD = DISPDD0 + rho*SOLN;
        DISPTD = DISPTD0 + rho*SOLN;           
    end
     
   
end