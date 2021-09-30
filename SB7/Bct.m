function [BI13,BI23]=Bct(XYZ,XI,I)
% Calcul les cisaillement transverse en un point
         [J1,J2,J3]=colomnejacobien(XYZ,XI);
         DN=dshapePrisme(XI);
         BI13=J1'*DN(3,I)+J3'*DN(1,I);
         BI23=J2'*DN(3,I)+J3'*DN(2,I);
end