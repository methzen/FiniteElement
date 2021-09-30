function [EcI13,EcI23]=Ect(DSP,XYZ,XI)
% Calcul les cisaillement transverse en un point
         [J1,J2,J3]=colomnejacobien(XYZ,XI);
         [D1,D2,D3]=colomneDefconv(DSP,XI);
         EcI13=dot(J1,D3)+dot(D1,J3)+dot(D1,D3);
         EcI23=dot(J3,D2)+dot(D3,J2)+dot(D3,D2);
end