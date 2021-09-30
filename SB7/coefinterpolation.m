function Bm=coefinterpolation(XYZ,XI,zeta)
r=XI(1);
s=XI(2);
r1=0.5-0.5/sqrt(3);
s1=r1;
r2=0.5+0.5/sqrt(3);
s2=r2;
%% rr
B1rr=zeros(1,18);
XI=[r1 0 zeta];
      [J1,J2,J3]=colomnejacobien(XYZ,XI);
      DN=dshapePrisme(XI);
      for I=1:6
         COL=3*(I-1)+1:3*(I-1)+3;
         B1rr(:,COL)=J1'*DN(1,I);
      end
      %%
B2rr=zeros(1,18);
XI=[r2 0 zeta];
      [J1,J2,J3]=colomnejacobien(XYZ,XI);
      DN=dshapePrisme(XI);
      for I=1:6
         COL=3*(I-1)+1:3*(I-1)+3;
         B2rr(:,COL)=J1'*DN(1,I);
      end
      %%
Bcrr=zeros(1,18);
XI=[r1 1/(sqrt(3)) zeta];
      [J1,J2,J3]=colomnejacobien(XYZ,XI);
      DN=dshapePrisme(XI);
      for I=1:6
         COL=3*(I-1)+1:3*(I-1)+3;
         Bcrr(:,COL)=J1'*DN(1,I);
      end
%% ss
B1ss=zeros(1,18);
XI=[0 s1 zeta];
      [J1,J2,J3]=colomnejacobien(XYZ,XI);
      DN=dshapePrisme(XI);
      for I=1:6
         COL=3*(I-1)+1:3*(I-1)+3;
         B1ss(:,COL)=J2'*DN(2,I);
      end      
%%      
B2ss=zeros(1,18);
XI=[0 s2 zeta];
      [J1,J2,J3]=colomnejacobien(XYZ,XI);
      DN=dshapePrisme(XI);
      for I=1:6
         COL=3*(I-1)+1:3*(I-1)+3;
         B2ss(:,COL)=J2'*DN(2,I);
      end
%%      
Bcss=zeros(1,18);
XI=[sqrt(1/3) s1 zeta];
      [J1,J2,J3]=colomnejacobien(XYZ,XI);
      DN=dshapePrisme(XI);
      for I=1:6
         COL=3*(I-1)+1:3*(I-1)+3;
         Bcss(:,COL)=J2'*DN(2,I);
      end
%% qq
B1qq=zeros(1,18);
XI=[r2 s1 zeta];
      [J1,J2,J3]=colomnejacobien(XYZ,XI);
      DN=dshapePrisme(XI);
      for I=1:6
         COL=3*(I-1)+1:3*(I-1)+3;
         B1qq(:,COL)=0.5*(J1'*DN(1,I)+J2'*DN(2,I))-(J1'*DN(2,I)+J2'*DN(1,I));
      end      
%%      
B2qq=zeros(1,18);
XI=[r1 s2 zeta];
      [J1,J2,J3]=colomnejacobien(XYZ,XI);
      DN=dshapePrisme(XI);
      for I=1:6
         COL=3*(I-1)+1:3*(I-1)+3;
         B2qq(:,COL)=0.5*(J1'*DN(1,I)+J2'*DN(2,I))-(J1'*DN(2,I)+J2'*DN(1,I));
      end
%%      
Bcqq=zeros(1,18);
XI=[r1 s1 zeta];
      [J1,J2,J3]=colomnejacobien(XYZ,XI);
      DN=dshapePrisme(XI);
      for I=1:6
         COL=3*(I-1)+1:3*(I-1)+3;
         Bcqq(:,COL)=0.5*(J1'*DN(1,I)+J2'*DN(2,I))-(J1'*DN(2,I)+J2'*DN(1,I));
      end
%%
a1=0.5*(B1rr+B2rr)-sqrt(3)*0.5*(B2rr-B1rr);
b1=sqrt(3)*(B2rr-B1rr);
c1=sqrt(3)*(Bcrr-a1-b1*r1);

a2=0.5*(B1ss+B2ss)-sqrt(3)*0.5*(B2ss-B1ss);
c2=sqrt(3)*(B2ss-B1ss);
b2=sqrt(3)*(Bcss-a2-c2*s1); 

a3=0.5*(B1qq+B2qq)-sqrt(3)*0.5*(B2qq-B1qq);
b3=-sqrt(3)*(B2qq-B1qq);
c3=sqrt(3)*(Bcqq-a3-b3*r1); 

Bm=zeros(3,18);
Bm(1:2,:)=[a1+r*b1+s*c1;
           a2+r*b2+s*c2];
    %a3+r*b3+c3*(1-r-s)
Bm(3,:)=0.5*(a1+r*b1+s*c1+a2+r*b2+s*c2)-(a3+r*b3+c3*(1-r-s));
end