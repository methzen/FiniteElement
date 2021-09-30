%%
%%calcul de la matrice jacobienne 
%%calcul�e � l'origine au noeud 1 du triangle central
%%J en ksi=0, eta=0, zeta=0
function [invJ]=MatricesJJm1Origine(xi,yi,zi,xj,yj,zj,xk,yk,zk...
          ,xl,yl,zl,xm,ym,zm,xn,yn,zn)
%%J en ksi=0, eta=0, zeta=0
      matJ(1,1)=(-xi+xj   -xl+xm   )/2.0;
      matJ(1,2)=(-yi+yj   -yl+ym   )/2.0;
      matJ(1,3)=(-zi+zj   -zl+zm   )/2.0;
      matJ(2,1)=(-xi   +xk-xl   +xn)/2.0;
      matJ(2,2)=(-yi   +yk-yl   +yn)/2.0;
      matJ(2,3)=(-zi   +zk-zl   +zn)/2.0;
      matJ(3,1)=(-xi      +xl      )/2.0;
      matJ(3,2)=(-yi      +yl      )/2.0;
      matJ(3,3)=(-zi      +zl      )/2.0;

%%determinant de J
      detJ=     (matJ(2,2)*matJ(3,3)-matJ(3,2)*matJ(2,3))*matJ(1,1);
      detJ=detJ-(matJ(1,2)*matJ(3,3)-matJ(3,2)*matJ(1,3))*matJ(2,1);
      detJ=detJ+(matJ(1,2)*matJ(2,3)-matJ(2,2)*matJ(1,3))*matJ(3,1);

%%matrice "J moins 1"
      invJ=inv(matJ);
   
end