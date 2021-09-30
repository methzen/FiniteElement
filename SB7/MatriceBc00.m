%%
%%pareil pour cisaillement transverse
%%mais consid�r� nul car la seule raideur pour le cisaillement transverse 
%%est introduite dans la stabilisation
function Bct=MatriceBc0(XYZ)

XYZ1(1,:)=XYZ(3,:);
XYZ1(2,:)=XYZ(1,:);
XYZ1(3,:)=XYZ(2,:);
XYZ1(4,:)=XYZ(6,:);
XYZ1(5,:)=XYZ(4,:);
XYZ1(6,:)=XYZ(5,:);
XYZ=XYZ1;

      xi=XYZ(1,1);yi=XYZ(1,2);zi=XYZ(1,3);
      xj=XYZ(2,1);yj=XYZ(2,2);zj=XYZ(2,3);
      xk=XYZ(3,1);yk=XYZ(3,2);zk=XYZ(3,3);
      xl=XYZ(4,1);yl=XYZ(4,2);zl=XYZ(4,3);
      xm=XYZ(5,1);ym=XYZ(5,2);zm=XYZ(5,3);
      xn=XYZ(6,1);yn=XYZ(6,2);zn=XYZ(6,3);
      %%
      h1=sqrt((xl-xi)*(xl-xi)+(yl-yi)*(yl-yi)+(zl-zi)*(zl-zi));
      h2=sqrt((xm-xj)*(xm-xj)+(ym-yj)*(ym-yj)+(zm-zj)*(zm-zj));
      h3=sqrt((xn-xk)*(xn-xk)+(yn-yk)*(yn-yk)+(zn-zk)*(zn-zk));

      [g304,g305,g306]=vecteurs_g30(h1,h2,h3,XYZ);
      
%%
           [V1xx,V1xy,V1xz,V1yx,V1yy,V1yz,V1zx,V1zy,V1zz,...
            V2xx,V2xy,V2xz,V2yx,V2yy,V2yz,V2zx,V2zy,V2zz,...
            V3xx,V3xy,V3xz,V3yx,V3yy,V3yz,V3zx,V3zy,V3zz]...
            =frame_orth_V0(XYZ,h1,h2,h3);

%       z1=(zi+zl)/2.0;
%       z2=(zj+zm)/2.0;
%       z3=(zk+zn)/2.0;      
%       g10(1)=x2-x1;g10(2)=y2-y1;g10(3)=z2-z1;
%       g20(1)=x3-x1;g20(2)=y3-y1;g20(3)=z3-z1;
%       g21(1)=x3-x2;g21(2)=y3-y2;g21(3)=z3-z2;
      [J1,J2,J3]=colomnejacobien(XYZ,[0 0 0]);
      g10=J1;
      g20=J2;
      g21=g20-g10;
%% calcul de matrices C
%     Ci=  Ci11  Ci12  Ci13
%          Ci21  Ci22  Ci23 

% initialisation
      C111=zeros(1,5);
      C112=zeros(1,5);
      C113=zeros(1,5);
      C121=zeros(1,5);
      C122=zeros(1,5);
      C123=zeros(1,5);

      C211=zeros(1,5);
      C212=zeros(1,5);
      C213=zeros(1,5);
      C221=zeros(1,5);
      C222=zeros(1,5);
      C223=zeros(1,5);
      
      C311=zeros(1,5);
      C312=zeros(1,5);
      C313=zeros(1,5);
      C321=zeros(1,5);
      C322=zeros(1,5);
      C323=zeros(1,5);
%      
      C111(1)=-1.0*g304(1); 
      C111(2)=-1.0*g304(2); 
      C111(3)=-1.0*g304(3);      
      C111(4)=-0.25*h1*(V1yx*g10(1)+V1yy*g10(2)+V1yz*g10(3));
      C111(5)=0.25*h1*(V1xx*g10(1)+V1xy*g10(2)+V1xz*g10(3));
      
      C112(1)=g304(1); 
      C112(2)=g304(2); 
      C112(3)=g304(3);     
      C112(4)=-0.25*h2*(V2yx*g10(1)+V2yy*g10(2)+V2yz*g10(3));
      C112(5)=0.25*h2*(V2xx*g10(1)+V2xy*g10(2)+V2xz*g10(3));

      C121(1)=-1.0*g306(1); 
      C121(2)=-1.0*g306(2); 
      C121(3)=-1.0*g306(3);      
      C121(4)=-0.25*h1*(V1yx*g10(1)+V1yy*g10(2)+V1yz*g10(3));
      C121(5)=+0.25*h1*(V1xx*g10(1)+V1xy*g10(2)+V1xz*g10(3));
      
      C123(1)=1.0*g306(1); 
      C123(2)=1.0*g306(2); 
      C123(3)=1.0*g306(3);      
      C123(4)=-0.25*h3*(V3yx*g20(1)+V3yy*g20(2)+V3yz*g20(3));
      C123(5)=+0.25*h3*(V3xx*g20(1)+V3xy*g20(2)+V3xz*g20(3));      
      C1=[C111 C112 C113;
          C121 C122 C123];      
%      
     %C211=0
      C212(1)=-1.0*g305(1); 
      C212(2)=-1.0*g305(2); 
      C212(3)=-1.0*g305(3);       
      C212(4)=-0.25*h2*(V2yx*g21(1)+V2yy*g21(2)+V2yz*g21(3));
      C212(5)=+0.25*h2*(V2xx*g21(1)+V2xy*g21(2)+V2xz*g21(3));

      C213(1)=g305(1); 
      C213(2)=g305(2); 
      C213(3)=g305(3);      
      C213(4)=-0.25*h3*(V3yx*g21(1)+V3yy*g21(2)+V3yz*g21(3));
      C213(5)=0.25*h3*(V3xx*g21(1)+V3xy*g21(2)+V3xz*g21(3));
      
      C221(1)=1.0*g304(1); 
      C221(2)=1.0*g304(2); 
      C221(3)=1.0*g304(3);       
      C221(4)=0.25*h1*(V1yx*g10(1)+V1yy*g10(2)+V1yz*g10(3));
      C221(5)=-0.25*h1*(V1xx*g10(1)+V1xy*g10(2)+V1xz*g10(3));

      C222(1)=-1.0*g304(1); 
      C222(2)=-1.0*g304(2); 
      C222(3)=-1.0*g304(3);       
      C222(4)=+0.25*h2*(V2yx*g10(1)+V2yy*g10(2)+V2yz*g10(3));
      C222(5)=-0.25*h2*(V2xx*g10(1)+V2xy*g10(2)+V2xz*g10(3));
     %C223=0
      C2=[C211 C212 C213;
          C221 C222 C223];      
      
      C311(1)=g306(1); 
      C311(2)=g306(2); 
      C311(3)=g306(3);      
      C311(4)=0.25*h1*(V1yx*g20(1)+V1yy*g20(2)+V1yz*g20(3));
      C311(5)=-0.25*h1*(V1xx*g20(1)+V1xy*g20(2)+V1xz*g20(3));
     %C312=0
     
      C313(1)=-1.0*g306(1); 
      C313(2)=-1.0*g306(2); 
      C313(3)=-1.0*g306(3);      
      C313(4)=+0.25*h3*(V3yx*g20(1)+V3yy*g20(2)+V3yz*g20(3));
      C313(5)=-0.25*h3*(V3xx*g20(1)+V3xy*g20(2)+V3xz*g20(3));
     %C321=0
      C322(1)=g305(1); 
      C322(2)=g305(2); 
      C322(3)=g305(3);      
      C322(4)=+0.25*h2*(V2yx*g21(1)+V2yy*g21(2)+V2yz*g21(3));
      C322(5)=-0.25*h2*(V2xx*g21(1)+V2xy*g21(2)+V2xz*g21(3));

      C323(1)=-g305(1); 
      C323(2)=-g305(2); 
      C323(3)=-g305(3);      
      C323(4)=+0.25*h3*(V3yx*g21(1)+V3yy*g21(2)+V3yz*g21(3));
      C323(5)=-0.25*h3*(V3xx*g21(1)+V3xy*g21(2)+V3xz*g21(3));
      C3=[C311 C312 C313;
          C321 C322 C323];      

%% assemblage de matrices C dans l'ordre de u1 v1 w1 u2 v2 w2 u3 
%% v3 w3 deta11 deta12 deta21 deta22 deta31 deta32       
      D(1,1)=1.0; 
      D(1,2)=0.0; 
      D(2,1)=0.0; 
      D(2,2)=1.0;
      aux1=D*C1;                      
      D(1,1)=0.0; 
      D(1,2)=-1.0; 
      D(2,1)=1.0; 
      D(2,2)=-1.0;     
      aux2=D*C2;
      D(1,1)=-1.0; 
      D(1,2)=1.0; 
      D(2,1)=-1.0; 
      D(2,2)=0.0;
      aux3=D*C3; 
       T=matrice_trans_rot(h1,h2,h3,...
            V1xx,V1xy,V1xz,V1yx,V1yy,V1yz,V1zx,V1zy,V1zz,...
            V2xx,V2xy,V2xz,V2yx,V2yy,V2yz,V2zx,V2zy,V2zz,...
            V3xx,V3xy,V3xz,V3yx,V3yy,V3yz,V3zx,V3zy,V3zz) ;
Bct=(1/3)*(aux1+aux2+aux3)*T;
B1=Bct(:,1:3);B2=Bct(:,4:6);B3=Bct(:,7:9);B4=Bct(:,10:12);B5=Bct(:,13:15);
B6=Bct(:,16:18);
Bct=[B2 B3 B1 B5 B6 B4];
 
end