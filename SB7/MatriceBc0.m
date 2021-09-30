%%
%%pareil pour cisaillement transverse
%%mais consid�r� nul car la seule raideur pour le cisaillement transverse 
%%est introduite dans la stabilisation
function [Bc0,Bc1,Bc2]=MatriceBc0(ni,nj,nk,xi,yi,zi,xj,yj,zj,xk,yk,zk...
          ,xl,yl,zl,xm,ym,zm,xn,yn,zn,...
          Xi,Yi,Zi,Xj,Yj,Zj,Xk,Yk,Zk,...
          Xl,Yl,Zl,Xm,Ym,Zm,Xn,Yn,Zn,...
          xX,xY,xZ,yX,yY,yZ,zX,zY,zZ,...
          x1,y1,x2,y2,x3,y3)

      h1=sqrt((xl-xi)*(xl-xi)+(yl-yi)*(yl-yi)+(zl-zi)*(zl-zi));
      h2=sqrt((xm-xj)*(xm-xj)+(ym-yj)*(ym-yj)+(zm-zj)*(zm-zj));
      h3=sqrt((xn-xk)*(xn-xk)+(yn-yk)*(yn-yk)+(zn-zk)*(zn-zk));
      for p=1:2
         for q=1:18
            Bc0(p,q)=0.0;
            Bc1(p,q)=0.0;
            Bc2(p,q)=0.0;
         end
      end  
%%J en ksi=0, eta=0, zeta=0
      [invJo]=MatricesJJm1Origine(xi,yi,zi,xj,yj,zj,xk,yk,zk...
          ,xl,yl,zl,xm,ym,zm,xn,yn,zn);
%%calcul les trois vecteurs g30 pour les trois mi-points de c�t�s
      [g304,g305,g306]=vecteurs_g30(h1,h2,h3,...
           xi,yi,zi,xj,yj,zj,xk,yk,zk,...
           xl,yl,zl,xm,ym,zm,xn,yn,zn);
%%calcul de frame orthogonal de rotation Vi10 Vi20 Vi30 pour chaque noeud
      [V1xx,V1xy,V1xz,V1yx,V1yy,V1yz,V1zx,V1zy,V1zz,...
            V2xx,V2xy,V2xz,V2yx,V2yy,V2yz,V2zx,V2zy,V2zz,...
            V3xx,V3xy,V3xz,V3yx,V3yy,V3yz,V3zx,V3zy,V3zz]...
            =frame_orth_V0(xi,yi,zi,xj,yj,zj,xk,yk,zk,...
            xl,yl,zl,xm,ym,zm,xn,yn,zn,...
            h1,h2,h3);
%       [V1xx,V1xy,V1xz,V1yx,V1yy,V1yz,V1zx,V1zy,V1zz,...
%             V2xx,V2xy,V2xz,V2yx,V2yy,V2yz,V2zx,V2zy,V2zz,...
%             V3xx,V3xy,V3xz,V3yx,V3yy,V3yz,V3zx,V3zy,V3zz]...
%             =frame_orth_V0(Xi,Yi,Zi,Xj,Yj,Zj,Xk,Yk,Zk,...
%             Xl,Yl,Zl,Xm,Ym,Zm,Xn,Yn,Zn,...
%             h1,h2,h3);
      z1=(zi+zl)/2.0;
      z2=(zj+zm)/2.0;
      z3=(zk+zn)/2.0;      
      g10(1)=x2-x1;g10(2)=y2-y1;g10(3)=z2-z1;
      g20(1)=x3-x1;g20(2)=y3-y1;g20(3)=z3-z1;
      g21(1)=x3-x2;g21(2)=y3-y2;g21(3)=z3-z2;
%%calcul de matrices C
%%     Ci=  Ci11  Ci12  Ci13
%%          Ci21  Ci22  Ci23          
      C111(1)=-1.0*g304(1); C111(2)=-1.0*g304(2); C111(3)=-1.0*g304(3);      
      C111(4)=-0.25*h1*(V1yx*g10(1)+V1yy*g10(2)+V1yz*g10(3));
      C111(5)=0.25*h1*(V1xx*g10(1)+V1xy*g10(2)+V1xz*g10(3));
      
      C112(1)=g304(1); C112(2)=g304(2); C112(3)=g304(3);     
      C112(4)=-0.25*h2*(V2yx*g10(1)+V2yy*g10(2)+V2yz*g10(3));
      C112(5)=0.25*h2*(V2xx*g10(1)+V2xy*g10(2)+V2xz*g10(3));
            
      C212(1)=-1.0*g305(1); C212(2)=-1.0*g305(2); C212(3)=-1.0*g305(3);       
      C212(4)=-0.25*h2*(V2yx*g21(1)+V2yy*g21(2)+V2yz*g21(3));
      C212(5)=0.25*h2*(V2xx*g21(1)+V2xy*g21(2)+V2xz*g21(3));
      
      C213(1)=g305(1); C213(2)=g305(2); C213(3)=g305(3);      
      C213(4)=-0.25*h3*(V3yx*g21(1)+V3yy*g21(2)+V3yz*g21(3));
      C213(5)=0.25*h3*(V3xx*g21(1)+V3xy*g21(2)+V3xz*g21(3));
      
      C311(1)=g306(1); C311(2)=g306(2); C311(3)=g306(3);      
      C311(4)=0.25*h1*(V1yx*g20(1)+V1yy*g20(2)+V1yz*g20(3));
      C311(5)=-0.25*h1*(V1xx*g20(1)+V1xy*g20(2)+V1xz*g20(3));
      
      C313(1)=-1.0*g306(1); C313(2)=-1.0*g306(2); C313(3)=-1.0*g306(3);      
      C313(4)=0.25*h3*(V3yx*g20(1)+V3yy*g20(2)+V3yz*g20(3));
      C313(5)=-0.25*h3*(V3xx*g20(1)+V3xy*g20(2)+V3xz*g20(3));
      
      C113(1)=0.0; C113(2)=0.0; C113(3)=0.0; C113(4)=0.0; C113(5)=0.0;
      C211(1)=0.0; C211(2)=0.0; C211(3)=0.0; C211(4)=0.0; C211(5)=0.0;
      C312(1)=0.0; C312(2)=0.0; C312(3)=0.0; C312(4)=0.0; C312(5)=0.0;
      
      for i=1:5
         C121(i)=-C311(i); C122(i)=-C312(i); C123(i)=-C313(i);  
         C221(i)=-C111(i); C222(i)=-C112(i); C223(i)=-C113(i);
         C321(i)=-C211(i); C322(i)=-C212(i); C323(i)=-C213(i);   
      end      
      for i=1:3
      C1(1,i)=C111(i); C1(1,i+3)=C112(i); C1(1,i+6)=C113(i);
      C1(2,i)=C121(i); C1(2,i+3)=C122(i); C1(2,i+6)=C123(i);
      
      C2(1,i)=C211(i); C2(1,i+3)=C212(i); C2(1,i+6)=C213(i);
      C2(2,i)=C221(i); C2(2,i+3)=C222(i); C2(2,i+6)=C223(i);
      
      C3(1,i)=C311(i); C3(1,i+3)=C312(i); C3(1,i+6)=C313(i);
      C3(2,i)=C321(i); C3(2,i+3)=C322(i); C3(2,i+6)=C323(i);
      end
%% assemblage de matrices C dans l'ordre de u1 v1 w1 u2 v2 w2 u3 
%% v3 w3 deta11 deta12 deta21 deta22 deta31 deta32     
      for i=4:5 
      if(ni>nj)
      C111(i)=-C111(i);C211(i)=-C211(i);C311(i)=-C311(i);
      C121(i)=-C121(i);C221(i)=-C221(i);C321(i)=-C321(i);
      end
      if(nj>nk)
      C112(i)=-C112(i);C212(i)=-C212(i);C312(i)=-C312(i);
      C122(i)=-C122(i);C222(i)=-C222(i);C322(i)=-C322(i);
      end
      if(nk>ni)
      C113(i)=-C113(i);C213(i)=-C213(i);C313(i)=-C313(i);
      C123(i)=-C123(i);C223(i)=-C223(i);C323(i)=-C323(i);
      end 
      C1(1,i+6)=C111(i); C1(1,i+8)=C112(i); C1(1,i+10)=C113(i);
      C1(2,i+6)=C121(i); C1(2,i+8)=C122(i); C1(2,i+10)=C123(i);
      
      C2(1,i+6)=C211(i); C2(1,i+8)=C212(i); C2(1,i+10)=C213(i);
      C2(2,i+6)=C221(i); C2(2,i+8)=C222(i); C2(2,i+10)=C223(i);
      
      C3(1,i+6)=C311(i); C3(1,i+8)=C312(i); C3(1,i+10)=C313(i);
      C3(2,i+6)=C321(i); C3(2,i+8)=C322(i); C3(2,i+10)=C323(i);
      end   
      D(1,1)=1.0; D(1,2)=0.0; D(2,1)=0.0; D(2,2)=1.0;
%%calcul matrice Bd = D*C     
      for p=1:2
         for q=1:15
            aux=0.0;
            for r=1:2
               aux=aux+D(p,r)*C1(r,q);
            end
            Bd1(p,q)=aux;
         end
      end
      D(1,1)=0.0; D(1,2)=-1.0; D(2,1)=1.0; D(2,2)=-1.0;     
      for p=1:2
         for q=1:15
            aux=0.0;
            for r=1:2
               aux=aux+D(p,r)*C2(r,q);
            end
            Bd2(p,q)=aux;
         end
      end
      D(1,1)=-1.0; D(1,2)=1.0; D(2,1)=-1.0; D(2,2)=0.0;      
      for p=1:2
         for q=1:15
            aux=0.0;
            for r=1:2
               aux=aux+D(p,r)*C3(r,q);
            end
            Bd3(p,q)=aux;
         end
      end 
%%calcul matrice Bdt � la base global 
%%calcul matrice transformation T
      for p=1:15
         for q=1:18
            T(p,q)=0.0;
         end
      end
      for p=1:3
      T(3*p-2,3*p-2)=0.5*xX; T(3*p-2,3*p-1)=0.5*xY;
      T(3*p-2,3*p)=0.5*xZ;
      T(3*p-1,3*p-2)=0.5*yX; T(3*p-1,3*p-1)=0.5*yY;
      T(3*p-1,3*p)=0.5*yZ;
      T(3*p,3*p-2)=0.5*zX; T(3*p,3*p-1)=0.5*zY;
      T(3*p,3*p)=0.5*zZ;
      
      T(3*p-2,3*p+7)=0.5*xX; T(3*p-2,3*p+8)=0.5*xY;
      T(3*p-2,3*p+9)=0.5*xZ;
      T(3*p-1,3*p+7)=0.5*yX; T(3*p-1,3*p+8)=0.5*yY;
      T(3*p-1,3*p+9)=0.5*yZ;
      T(3*p,3*p+7)=0.5*zX; T(3*p,3*p+8)=0.5*zY;
      T(3*p,3*p+9)=0.5*zZ;
      end
      [Tr]=matrice_trans_rot(ni,nj,nk,h1,h2,h3,...
            V1xx,V1xy,V1xz,V1yx,V1yy,V1yz,V1zx,V1zy,V1zz,...
            V2xx,V2xy,V2xz,V2yx,V2yy,V2yz,V2zx,V2zy,V2zz,...
            V3xx,V3xy,V3xz,V3yx,V3yy,V3yz,V3zx,V3zy,V3zz,...
            xX,xY,xZ,yX,yY,yZ,zX,zY,zZ); 
%         Tr=zeros(6,18);  
      for p=1:6
         for q=1:18
            T(p+9,q)=Tr(p,q);
         end
      end
      
      for p=1:2
         for q=1:18
            aux=0.0;
            for r=1:15
               aux=aux+Bd1(p,r)*T(r,q);
            end
            Bdt1(p,q)=aux;
         end
      end
      for p=1:2
         for q=1:18
            aux=0.0;
            for r=1:15
               aux=aux+Bd2(p,r)*T(r,q);
            end
            Bdt2(p,q)=aux;
         end
      end  
      for p=1:2
         for q=1:18
            aux=0.0;
            for r=1:15
               aux=aux+Bd3(p,r)*T(r,q);
            end
            Bdt3(p,q)=aux;
         end
      end
      for p=1:2
         for q=1:18
            aux=0.0;
            for r=1:2
               aux=aux+invJo(p,r)*Bdt1(r,q);
            end
            Bc0(p,q)=aux;
         end
      end
      for p=1:2
         for q=1:18
            aux=0.0;
            for r=1:2
               aux=aux+invJo(p,r)*(Bdt3(r,q)-Bdt1(r,q));
            end
            Bc1(p,q)=aux;
         end
      end
      for p=1:2
         for q=1:18
            aux=0.0;
            for r=1:2
               aux=aux+invJo(p,r)*(Bdt2(r,q)-Bdt1(r,q));
            end
            Bc2(p,q)=aux;
         end
      end
end