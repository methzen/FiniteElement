function Ec=deformatEc(DSP,XYZ,XI)

% ICI LES COORDONNÉES SONT LES COORDONNÉES INITIALES
%
[J1,J2,J3]=colomnejacobien(XYZ,XI);
[D1,D2,D3]=colomneDefconv(DSP,XI);
Ec=zeros(6,1);
Ec(1)=dot(J1,D1)+0.5*dot(D1,D1);
Ec(2)=dot(J2,D2)+0.5*dot(D2,D2);
Ec(3)=dot(J3,D3)+0.5*dot(D3,D3);
Ec(4)=dot(J1,D2)+dot(D1,J2)+dot(D1,D2);
Ec(5)=dot(J2,D3)+dot(D2,J3)+dot(D2,D3);
Ec(6)=dot(J3,D1)+dot(D3,J1)+dot(D3,D1);

end