function Bm=Bmembrane(XYZ,XI)

Bm1=coefinterpolation(XYZ,XI,-1);
Bm2=coefinterpolation(XYZ,XI,+1);
Bm=0.5*(1-XI(3))*Bm1+0.5*(1+XI(3))*Bm2;
end