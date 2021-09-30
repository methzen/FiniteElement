function BN=Bpinch(XYZ,XI)

BNA=zeros(1,18);
      [J1,J2,J3A]=colomnejacobien(XYZ,[0 0 0]);
       DNA=dshapePrisme([0 0 0]);
BNB=zeros(1,18);
      [J1,J2,J3B]=colomnejacobien(XYZ,[1 0 0]);
       DNB=dshapePrisme([1 0 0]);
BNC=zeros(1,18);
      [J1,J2,J3C]=colomnejacobien(XYZ,[0 1 0]);
       DNC=dshapePrisme([0 1 0]);
      for I=1:6
         COL=3*(I-1)+1:3*(I-1)+3;
         BNA(:,COL)=J3A'*DNA(3,I);
         BNB(:,COL)=J3B'*DNB(3,I);
         BNC(:,COL)=J3C'*DNC(3,I);
      end
BN=(BNA+BNB+BNC)/3;
end