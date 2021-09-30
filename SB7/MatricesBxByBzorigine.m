%%
%%en ksi = 0 et eta = 0 ; "ksi" car xi existe
%%ligne 1 : <N,ksi> ; ligne 2 : <N,eta> ; ligne 3 : <N,zeta>
%%ligne 1:<N,x>; ligne 2:<N,y>; ligne 3:<N,z> - en xi=0, eta=0, zeta=0 -
function [Bx,By,Bz]=MatricesBxByBzorigine(invJ)

%%en ksi = 0 et eta = 0 ; "ksi" car xi existe
   zeta=0.0;
   b=(1.0-zeta)/2.0; t=(1.0+zeta)/2.0; u=1.0/2.0;

%%ligne 1 : <N,ksi> ; ligne 2 : <N,eta> ; ligne 3 : <N,zeta>
   Bksi(1,1)=-b;  Bksi(2,1)=-b;  Bksi(3,1)=-u;
   Bksi(1,2)= b;  Bksi(2,2)=0.0; Bksi(3,2)=0.0;
   Bksi(1,3)=0.0; Bksi(2,3)= b;  Bksi(3,3)=0.0;
   Bksi(1,4)=-t;  Bksi(2,4)=-t;  Bksi(3,4)= u;
   Bksi(1,5)= t;  Bksi(2,5)=0.0; Bksi(3,5)=0.0;
   Bksi(1,6)=0.0; Bksi(2,6)= t;  Bksi(3,6)=0.0;

%%ligne 1:<N,x>; ligne 2:<N,y>; ligne 3:<N,z> - en xi=1/3, eta=1/3, zeta=0 -
   for j=1:6
         Bx(j)=invJ(1,1)*Bksi(1,j)+invJ(1,2)*Bksi(2,j)+invJ(1,3)*Bksi(3,j);
         By(j)=invJ(2,1)*Bksi(1,j)+invJ(2,2)*Bksi(2,j)+invJ(2,3)*Bksi(3,j);
         Bz(j)=invJ(3,1)*Bksi(1,j)+invJ(3,2)*Bksi(2,j)+invJ(3,3)*Bksi(3,j);
   end
end