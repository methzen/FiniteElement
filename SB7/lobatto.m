    
% choix du nompbre de poits de gauss
function [ZG,WG]=lobatto(npg)

switch npg
  case -3    %{debut points de GAUSS}
      ZG=zeros(1,3);
      WG=zeros(1,3);      
        ZG(3)=sqrt(3.0/5.0);
    	ZG(2)=0.0;
        ZG(1)=-ZG(3);
        WG(3)=5.0/9.0;
    	WG(2)=8.0/9.0;
        WG(1)=WG(3);
    case 3
     ZG=[-1.0 0.0 1.0];
     WG=[1.0/3.0 4.0/3.0 1.0/3.0];
    case 5
     ZG=[-1 -sqrt(3/7) 0 sqrt(3/7)  1];
     WG=[0.1 49.0/90.0 32.0/45.0 49.0/90.0 0.1];
    case 7   
    ZG=[-1.0 -0.83022390 -0.46884879 0.0 0.46884879 0.83022390 1.0];
    WG=[0.04761904 0.2768260 0.43174538 0.48761904 0.43174538 0.27682604 0.04761904] ;  
end