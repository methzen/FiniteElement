%%
%%coordonn�es dans le rep�re local des six noeuds 
%%+ 1 2 3 les trois noeuds du plan central
%%ca donne l'aire de l'�l�ment central et l'�paisseur � 
%%consid�rer pour le prisme
function [xi,yi,zi,xj,yj,zj,xk,yk,zk,xl,yl,zl,xm,ym,zm,xn,yn,zn,...
         x1,y1,x2,y2,x3,y3,aire2D,h]...
         =CoordonneesLocalesij(Xi,Yi,Zi,Xj,Yj,Zj,Xk,Yk,Zk,...
         Xl,Yl,Zl,Xm,Ym,Zm,Xn,Yn,Zn,xX,xY,xZ,yX,yY,yZ,zX,zY,zZ)

      xi=xX*Xi+xY*Yi+xZ*Zi;
      yi=yX*Xi+yY*Yi+yZ*Zi;
      zi=zX*Xi+zY*Yi+zZ*Zi;
      xj=xX*Xj+xY*Yj+xZ*Zj;
      yj=yX*Xj+yY*Yj+yZ*Zj;
      zj=zX*Xj+zY*Yj+zZ*Zj;
      xk=xX*Xk+xY*Yk+xZ*Zk;
      yk=yX*Xk+yY*Yk+yZ*Zk;
      zk=zX*Xk+zY*Yk+zZ*Zk;

      xl=xX*Xl+xY*Yl+xZ*Zl;
      yl=yX*Xl+yY*Yl+yZ*Zl;
      zl=zX*Xl+zY*Yl+zZ*Zl;
      xm=xX*Xm+xY*Ym+xZ*Zm;
      ym=yX*Xm+yY*Ym+yZ*Zm;
      zm=zX*Xm+zY*Ym+zZ*Zm;
      xn=xX*Xn+xY*Yn+xZ*Zn;
      yn=yX*Xn+yY*Yn+yZ*Zn;
      zn=zX*Xn+zY*Yn+zZ*Zn;

      x1=(xi+xl)/2.0; x2=(xj+xm)/2.0; x3=(xk+xn)/2.0;
      y1=(yi+yl)/2.0; y2=(yj+ym)/2.0; y3=(yk+yn)/2.0;
      x3=x3-x1; x2=x2-x1; x1=0.0;
      y3=y3-y1; y2=0.0;   y1=0.0;
      aire2D=abs(x2*y3-x3*y2)/2.0;
      h=abs(zl+zm+zn-zi-zj-zk)/3.0;

end 