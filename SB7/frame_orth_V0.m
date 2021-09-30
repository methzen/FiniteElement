function [V1xx,V1xy,V1xz,V1yx,V1yy,V1yz,V1zx,V1zy,V1zz,...
            V2xx,V2xy,V2xz,V2yx,V2yy,V2yz,V2zx,V2zy,V2zz,...
            V3xx,V3xy,V3xz,V3yx,V3yy,V3yz,V3zx,V3zy,V3zz]...
            =frame_orth_V0(XYZ,h1,h2,h3)
%%        
      xi=XYZ(1,1);yi=XYZ(1,2);zi=XYZ(1,3);
      xj=XYZ(2,1);yj=XYZ(2,2);zj=XYZ(2,3);
      xk=XYZ(3,1);yk=XYZ(3,2);zk=XYZ(3,3);
      xl=XYZ(4,1);yl=XYZ(4,2);zl=XYZ(4,3);
      xm=XYZ(5,1);ym=XYZ(5,2);zm=XYZ(5,3);
      xn=XYZ(6,1);yn=XYZ(6,2);zn=XYZ(6,3);
       
      x1=(xi+xl)/2.0; x2=(xj+xm)/2.0; x3=(xk+xn)/2.0;
      y1=(yi+yl)/2.0; y2=(yj+ym)/2.0; y3=(yk+yn)/2.0;  
      z1=(zi+zl)/2.0;
      z2=(zj+zm)/2.0;
      z3=(zk+zn)/2.0;
%%v10           
      V1zx=(xl-xi)/h1;
      V1zy=(yl-yi)/h1;
      V1zz=(zl-zi)/h1;

      V1yx=V1zy*(z2-z1)-V1zz*(y2-y1);
      V1yy=V1zz*(x2-x1)-V1zx*(z2-z1);
      V1yz=V1zx*(y2-y1)-V1zy*(x2-x1);
      Ny=sqrt(V1yx*V1yx+V1yy*V1yy+V1yz*V1yz);
      V1yx=V1yx/Ny;
      V1yy=V1yy/Ny;
      V1yz=V1yz/Ny;

      V1xx=V1yy*V1zz-V1zy*V1yz;
      V1xy=V1yz*V1zx-V1zz*V1yx;
      V1xz=V1yx*V1zy-V1zx*V1yy;
      Nx=sqrt(V1xx*V1xx+V1xy*V1xy+V1xz*V1xz);
      V1xx=V1xx/Nx;
      V1xy=V1xy/Nx;
      V1xz=V1xz/Nx;
%%v20           
      V2zx=(xm-xj)/h2;
      V2zy=(ym-yj)/h2;
      V2zz=(zm-zj)/h2;

      V2yx=V2zy*(z3-z2)-V2zz*(y3-y2);
      V2yy=V2zz*(x3-x2)-V2zx*(z3-z2);
      V2yz=V2zx*(y3-y2)-V2zy*(x3-x2);
      Ny=sqrt(V2yx*V2yx+V2yy*V2yy+V2yz*V2yz);
      V2yx=V2yx/Ny;
      V2yy=V2yy/Ny;
      V2yz=V2yz/Ny;

      V2xx=V2yy*V2zz-V2zy*V2yz;
      V2xy=V2yz*V2zx-V2zz*V2yx;
      V2xz=V2yx*V2zy-V2zx*V2yy;
      Nx=sqrt(V2xx*V2xx+V2xy*V2xy+V2xz*V2xz);
      V2xx=V2xx/Nx;
      V2xy=V2xy/Nx;
      V2xz=V2xz/Nx;
%%v30           
      V3zx=(xn-xk)/h3;
      V3zy=(yn-yk)/h3;
      V3zz=(zn-zk)/h3;

      V3yx=V3zy*(z1-z3)-V3zz*(y1-y3);
      V3yy=V3zz*(x1-x3)-V3zx*(z1-z3);
      V3yz=V3zx*(y1-y3)-V3zy*(x1-x3);
      Ny=sqrt(V3yx*V3yx+V3yy*V3yy+V3yz*V3yz);
      V3yx=V3yx/Ny;
      V3yy=V3yy/Ny;
      V3yz=V3yz/Ny;

      V3xx=V3yy*V3zz-V3zy*V3yz;
      V3xy=V3yz*V3zx-V3zz*V3yx;
      V3xz=V3yx*V3zy-V3zx*V3yy;
      Nx=sqrt(V3xx*V3xx+V3xy*V3xy+V3xz*V3xz);
      V3xx=V3xx/Nx;
      V3xy=V3xy/Nx;
      V3xz=V3xz/Nx;              
               
end 