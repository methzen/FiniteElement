function linesearch(PROP,SOLN,FEXT,FREEDOF,XYZ,Connec,rtu0)
global DISPTD DISPDD

% on enregistre les déplacements precedente
    DISPTD0=DISPTD;
    DISPDD0=DISPDD;
    msearch=5;
    searc=0.1;
    DISPTD = DISPTD + 1*SOLN;
    RESIDU=forceinterne(PROP,DISPTD,FEXT,XYZ,Connec);
    rtu=full(dot(SOLN(FREEDOF),RESIDU(FREEDOF)));
    
    if abs(rtu) < abs(rtu0*searc)
        DISPDD = DISPDD0 + SOLN;
        DISPTD = DISPTD0 + SOLN;    
    else
          rho0  = 0;
          rho   = 1;   
          nsear = 0;  
          while((abs(rtu) >= abs(rtu0*searc) )&& (nsear < msearch))
              
                  [~,rho, ~] = search(rho0,rho,rtu0,rtu);  
                  DISPTD = DISPTD0 + rho*SOLN;
                  RESIDU=forceinterne(PROP,DISPTD,FEXT, XYZ,Connec);
                  rtu=full(dot(SOLN(FREEDOF),RESIDU(FREEDOF)));
                

                  nsear = nsear + 1 ;
          end 
        DISPDD = DISPDD0 + rho*SOLN;
        DISPTD = DISPTD0 + rho*SOLN;           
    end
     
   
end



