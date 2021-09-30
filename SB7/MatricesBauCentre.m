%%
%%calcul_Bg: bool�en pour savoir si oui ou non on calcul les matrices b g�om�triques
%%Xi,Yi,Zi,Xj,Yj,Zj,Xk,Yk,Zk   
%%Xl,Yl,Zl,Xm,Ym,Zm,Xn,Yn,Zn
%%--> coordonn�es dans le rep�re global des 6 noeuds du prismes
%%calcul au barycentre
%%invJo
%%Bm0,Bb0,Bgu,Bgv,Bgw --> Bb0 c'est celui du DKT6
%%Bc0 --> nul
%%Bp0--> pour le pincement 
%%Bsx,Bsy,Bsz --> pour la stabilisation
function [aire2D,h,volume,invJo,Bm0,Bb0,Bgu,Bgv,Bgw,...
            Bc0,Bc1,Bc2,Bp0,Bsx,Bsy,Bsz]=MatricesBauCentre(ni,nj,nk,XYZ)
       Xi=XYZ(1,1);
       Yi=XYZ(1,2);
       Zi=XYZ(1,3);
       Xj=XYZ(2,1);
       Yj=XYZ(2,2);
       Zj=XYZ(2,3);
       
       Xk=XYZ(3,1);
       Yk=XYZ(3,2);
       Zk=XYZ(3,3);
       
        Xl=XYZ(4,1);
        Yl=XYZ(4,2);
        Zl=XYZ(4,3);
        
        Xm=XYZ(5,1);
        Ym=XYZ(5,2);
        Zm=XYZ(5,3);
        
        Xn=XYZ(6,1);
        Yn=XYZ(6,2);
        Zn=XYZ(6,3);
%% calcul du rep�re local dans le plan de l'�l�ment central avec i pour origine
        [xX,xY,xZ,yX,yY,yZ,zX,zY,zZ]...
         =AxesLocauxij(Xi,Yi,Zi,Xj,Yj,Zj,Xk,Yk,Zk,Xl,Yl,Zl,Xm,Ym,Zm,Xn,Yn,Zn);
 
%%
%%coordonn�es dans le rep�re local des six noeuds 
%%+ 1 2 3 les trois noeuds du plan central
%%ca donne l'aire de l'�l�ment central et l'�paisseur � 
%%consid�rer pour le prisme
        [xi,yi,zi,xj,yj,zj,xk,yk,zk,xl,yl,zl,xm,ym,zm,xn,yn,zn...
            ,x1,y1,x2,y2,x3,y3,aire2D,h]...
         =CoordonneesLocalesij(Xi,Yi,Zi,Xj,Yj,Zj,Xk,Yk,Zk...
         ,Xl,Yl,Zl,Xm,Ym,Zm,Xn,Yn,Zn,xX,xY,xZ,yX,yY,yZ,zX,zY,zZ);
%% 
volume=aire2D*h;
%%
%%debut flexion
%%
%%appel des deux proc�dures pr�c�dentes Bt_Bw
[Bt,Bw]=matricesBtBw(ni,nj,nk,x1,y1,x2,y2,x3,y3,aire2D);
%%
%%donne les trois rotations aux milieux des cot�s du DKT6 en fonction des 
%%translations dans le rep�re global des 6 noeuds
[Tt]=TransformerTransRot(ni,nj,nk,Xi,Yi,Zi,Xj,Yj,Zj,Xk,Yk,Zk,...
             Xl,Yl,Zl,Xm,Ym,Zm,Xn,Yn,Zn);
         
         
%%
Bbt=Bt*Tt;
%  for p=1:3
%       for q=1:18
%           aux=0.0;
%           for r=1:3
%               aux=aux+Bt(p,r)*Tt(r,q);
%           Bbt(p,q)=aux;
%           end
%       end
%  end
  
  for p=1:3
   Bbw(p,1)=0.5*Bw(p,1)*zX; Bbw(p,2)=0.5*Bw(p,1)*zY; Bbw(p,3)=0.5*Bw(p,1)*zZ;
   Bbw(p,10)=Bbw(p,1);      Bbw(p,11)=Bbw(p,2);      Bbw(p,12)=Bbw(p,3);
   Bbw(p,4)=0.5*Bw(p,2)*zX; Bbw(p,5)=0.5*Bw(p,2)*zY; Bbw(p,6)=0.5*Bw(p,2)*zZ;
   Bbw(p,13)=Bbw(p,4);      Bbw(p,14)=Bbw(p,5);      Bbw(p,15)=Bbw(p,6);
   Bbw(p,7)=0.5*Bw(p,3)*zX; Bbw(p,8)=0.5*Bw(p,3)*zY; Bbw(p,9)=0.5*Bw(p,3)*zZ;
   Bbw(p,16)=Bbw(p,7);      Bbw(p,17)=Bbw(p,8);      Bbw(p,18)=Bbw(p,9);
  end
  
  Bb0=Bbt+Bbw;
 % for p=1:3
  %   for q=1:18
   %      Bb0(p,q)=Bbt(p,q)+Bbw(p,q);
    % end
 %end
%%fin flexion
%%
%%J en ksi=1/3, eta=1/3, zeta=0
%%matrice "J moins 1"
[detJ,matJ,invJ]=MatricesJJm1auCentre(xi,yi,zi,xj,yj,zj,xk,yk,zk...
          ,xl,yl,zl,xm,ym,zm,xn,yn,zn);
%%
[bx,by,bz]=MatricesBxByBzauCentre(invJ);
%%
%%B membrane
%%en input ce sont les translations dans le rep�re global et renvoit l'accroissement 
%%de la d�formation membranaire 
%%n'est pas utilis� pour la stabilisation
[Bm0]=MatriceBm0(xX,xY,xZ,yX,yY,yZ,bx,by);
%%
%%parail pour le pincement
[Bp0]=MatriceBp0(zX,zY,zZ,bz);
%%matrice Bp zeta
%[Bpzeta]=MatriceBpzeta(zX,zY,zZ,h,5);

%%
%%pareil pour cisaillement transverse
%%mais consid�r� nul car la seule raideur pour le cisaillement transverse 
%%est introduite dans la stabilisation
[Bc0,Bc1,Bc2]=MatriceBc0(ni,nj,nk,xi,yi,zi,xj,yj,zj,xk,yk,zk...
          ,xl,yl,zl,xm,ym,zm,xn,yn,zn,...
          Xi,Yi,Zi,Xj,Yj,Zj,Xk,Yk,Zk,...
          Xl,Yl,Zl,Xm,Ym,Zm,Xn,Yn,Zn,...
          xX,xY,xZ,yX,yY,yZ,zX,zY,zZ,...
          x1,y1,x2,y2,x3,y3);
%            C=Bc0'*Bc0;
%           M=Bm0'*Bm0;
%           Bm0=zeros(3,18);Bp0=zeros(1,18);Bb0=zeros(3,18);
%         Bc0=zeros(2,18);Bc1=zeros(2,18);Bc2=zeros(2,18);
%%
%%matrice B g�om�trique 
%%en rep�re local
Bgu(3,18)=0;
Bgv(3,18)=0;
Bgw(3,18)=0;
%  if (calculBg~=0)
%  [Bgu,Bgv,Bgw]=MatricesBg(xX,xY,xZ,yX,yY,yZ,zX,zY,zZ,bx,by,bz);
%  end
%%
%%stabilisation
%%
%%calcul de la matrice jacobienne 
%%calcul�e � l'origine au noeud 1 du triangle central
%%J en ksi=0, eta=0, zeta=0
[invJo]=MatricesJJm1Origine(xi,yi,zi,xj,yj,zj,xk,yk,zk...
          ,xl,yl,zl,xm,ym,zm,xn,yn,zn);
%%
%%en ksi = 0 et eta = 0 ; "ksi" car xi existe
[Bx,By,Bz]=MatricesBxByBzorigine(invJo);
%%
%%calcul des vecteurs gamma de la stabilisation
%%il y a normalement 4 pour l'hexa�dre mais l� on n'en a que deux pour un prisme
[Vgamma]=VecteursGamma(xi,yi,zi,xj,yj,zj,xk,yk,zk...
          ,xl,yl,zl,xm,ym,zm,xn,yn,zn,Bx,By,Bz);
%%
%%matrice gradient pour la stabilisation
%%Bsx pour q1x,q1y; Bsy pour q1y,q2y; Bsz pour q1z,q2z "d�formations" de stabilisation
[Bsx,Bsy,Bsz]=MatricesBs(xX,xY,xZ,yX,yY,yZ,zX,zY,zZ,Vgamma);


end