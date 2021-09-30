%%
%%calcul du rep�re d'une face lat�riales en son centre
%%-Xi,Yi,Zi,Xj,Yj,Zj,Xm,Ym,Zm,Xl,Yl,Zl  
%%--> coordonn�es des noeuds dans le rep�re global
%%-xXf,xYf,xZf,yXf,yYf,yZf,zXf,zYf,zZf
%%--> cosinus directeur 
%%-xif,yif,xjf,yjf,xlf,ylf,xmf,ymf,dAf
%--> coordonn�es dans ce rep�re local des 4 neuds de la face lat�rale 
function [xXf,xYf,xZf,yXf,yYf,yZf,zXf,zYf,zZf...
          xif,yif,xjf,yjf,xlf,ylf,xmf,ymf,dAf...
          ]=RepereFaceMilieuxij...
          (Xi,Yi,Zi,Xj,Yj,Zj,Xm,Ym,Zm,Xl,Yl,Zl)
      Xab=(Xj+Xm)/2.0-(Xi+Xl)/2.0;
      Yab=(Yj+Ym)/2.0-(Yi+Yl)/2.0;
      Zab=(Zj+Zm)/2.0-(Zi+Zl)/2.0;

      Xcd=(Xl+Xm)/2.0-(Xi+Xj)/2.0;
      Ycd=(Yl+Ym)/2.0-(Yi+Yj)/2.0;
      Zcd=(Zl+Zm)/2.0-(Zi+Zj)/2.0;

      zXf=Yab*Zcd-Zab*Ycd;
      zYf=Zab*Xcd-Xab*Zcd;
      zZf=Xab*Ycd-Yab*Xcd;
      Nz=sqrt(zXf*zXf+zYf*zYf+zZf*zZf);
      zXf=zXf/Nz;
      zYf=zYf/Nz;
      zZf=zZf/Nz;

      Nx=sqrt(Xab*Xab+Yab*Yab+Zab*Zab);
      xXf=Xab/Nx;
      xYf=Yab/Nx;
      xZf=Zab/Nx;

      yXf=zYf*xZf-zZf*xYf;
      yYf=zZf*xXf-zXf*xZf;
      yZf=zXf*xYf-zYf*xXf;
      Ny=sqrt(yXf*yXf+yYf*yYf+yZf*yZf);
      yXf=yXf/Ny;
      yYf=yYf/Ny;
      yZf=yZf/Ny;

      xif=xXf*Xi+xYf*Yi+xZf*Zi;
      yif=yXf*Xi+yYf*Yi+yZf*Zi;
      xjf=xXf*Xj+xYf*Yj+xZf*Zj;
      yjf=yXf*Xj+yYf*Yj+yZf*Zj;
      xmf=xXf*Xm+xYf*Ym+xZf*Zm;
      ymf=yXf*Xm+yYf*Ym+yZf*Zm;
      xlf=xXf*Xl+xYf*Yl+xZf*Zl;
      ylf=yXf*Xl+yYf*Yl+yZf*Zl;

      dAf=abs(((xmf-xif)*(ylf-yjf)+(ymf-yif)*(xjf-xlf)))/2.0;    
end