function BN=calculB(XYZ,XI)

BN=zeros(6,18);
      [J1,J2,J3]=colomnejacobien(XYZ,XI);
      DN=dshapePrisme(XI);
      for I=1:6
         COL=3*(I-1)+1:3*(I-1)+3;
         BN(:,COL)=[J1'*DN(1,I);
                    J2'*DN(2,I);
                    J3'*DN(3,I);
                    J1'*DN(2,I)+J2'*DN(1,I);
                    J2'*DN(3,I)+J3'*DN(2,I);
                    J1'*DN(3,I)+J3'*DN(1,I)];
      end
end