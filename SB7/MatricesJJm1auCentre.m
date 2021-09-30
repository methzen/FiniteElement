%%
%%J en ksi=1/3, eta=1/3, zeta=0
%%matrice "J moins 1"
function [detJ,matJ,invJ]=MatricesJJm1auCentre(xi,yi,zi,xj,yj,zj,xk,yk,zk...
          ,xl,yl,zl,xm,ym,zm,xn,yn,zn)
%%J en ksi=1/3, eta=1/3, zeta=0
      matJ(1,1)=(-xi+xj   -xl+xm   )/2.0;
      matJ(1,2)=(-yi+yj   -yl+ym   )/2.0;
      matJ(1,3)=(-zi+zj   -zl+zm   )/2.0;
      matJ(2,1)=(-xi   +xk-xl   +xn)/2.0;
      matJ(2,2)=(-yi   +yk-yl   +yn)/2.0;
      matJ(2,3)=(-zi   +zk-zl   +zn)/2.0;
      matJ(3,1)=(-xi-xj-xk+xl+xm+xn)/6.0;
      matJ(3,2)=(-yi-yj-yk+yl+ym+yn)/6.0;
      matJ(3,3)=(-zi-zj-zk+zl+zm+zn)/6.0;

%%determinant de J
      detJ=     (matJ(2,2)*matJ(3,3)-matJ(3,2)*matJ(2,3))*matJ(1,1);
      detJ=detJ-(matJ(1,2)*matJ(3,3)-matJ(3,2)*matJ(1,3))*matJ(2,1);
      detJ=detJ+(matJ(1,2)*matJ(2,3)-matJ(2,2)*matJ(1,3))*matJ(3,1);

%%matrice "J moins 1"
%       invJ=Mat3InverseFast(matJ,detJ);
  invJ=inv(matJ);
  return  
end