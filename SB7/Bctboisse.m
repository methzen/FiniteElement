function Bct=Bctboisse(XYZ,zeta)
% Calcul les cisaillement transverse en un point
%%
% def noeud 1, m=[0.5 0.5 0]
B_n1=zeros(2,9);
XIm=[0.5 0.5 zeta];
         [J1,J2,J3]=colomnejacobien(XYZ,[0 0 zeta]);
         DN=dshapePrisme(XIm);
         g1=-J1+J2;
         g2=-J2;
         g3=J3;
         B_n1(1,1:3)=-g3'+DN(3,1)*g1';
         B_n1(1,4:6)=+g3'+DN(3,2)*g1';
         B_n1(1,7:9)=+0  +DN(3,3)*g1';
         B_n1(2,1:3)=-g3'+DN(3,1)*g2';
         B_n1(2,4:6)=+ 0 +DN(3,2)*g2';
         B_n1(2,7:9)=+g3'+DN(3,3)*g2'; 
D1=[0 -1;
    1 -1];
B_n1=D1*B_n1;
%%         
% def noeud 2, m=[0 0.5 0]
B_n2=zeros(2,9);
XIm=[0.0 0.5 zeta];
         [J1,J2,J3]=colomnejacobien(XYZ,[0 0 zeta]);
         DN=dshapePrisme(XIm);
         g1=-J2;
         g2=J1-J2;
         g3=J3;
         B_n2(1,1:3)=+ 0'+DN(3,1)*g1';
         B_n2(1,4:6)=-g3'+DN(3,2)*g1';
         B_n2(1,7:9)=+g3'+DN(3,3)*g1';
         
         B_n2(2,1:3)=+g3'+DN(3,1)*g2';
         B_n2(2,4:6)=-g3'+DN(3,2)*g2';
         B_n2(2,7:9)=+ 0 +DN(3,3)*g2';           
D2=[-1 1;
    -1 0];
B_n2=D2*B_n2;         
%%         
% def noeud 3, m=[0.5 0. 0]
B_n3=zeros(2,9);
XIm=[0.5 0. zeta];
         [J1,J2,J3]=colomnejacobien(XYZ,[0 0 zeta]);
         DN=dshapePrisme(XIm);
         g1=J1;
         g2=J2;
         g3=J3;
         B_n3(1,1:3)=+g3'+DN(3,1)*g1';
         B_n3(1,4:6)=+0'+DN(3,2)*g1';
         B_n3(1,7:9)=-g3'+DN(3,3)*g1';
         
         B_n3(2,1:3)=+g3'+DN(3,1)*g2';
         B_n3(2,4:6)=+0' +DN(3,2)*g2';
         B_n3(2,7:9)=-g3'+DN(3,3)*g2';
D3=[1 0;
    0 1];
B_n3=D3*B_n3;
%%
Bct=(B_n1+B_n2+B_n3)/3.0;
if zeta==1
    Bct=[zeros(2,9) Bct];
elseif zeta==-1
    Bct=[Bct zeros(2,9)];
end
end
