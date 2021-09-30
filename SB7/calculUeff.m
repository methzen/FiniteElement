function Ueff=calculUeff(S,E)

if E==0
    Ueff=0.0;
else
IDEN=[1 1 1 0 0 0]';
Sdev=S-(S(1)+S(2)+S(3))*IDEN/3;
Edev=E-(E(1)+E(2)+E(3))*IDEN/3;

Sdev3=[Sdev(1) Sdev(4) Sdev(6);
       Sdev(4) Sdev(2) Sdev(5);
       Sdev(6) Sdev(5) Sdev(3)];
Edev3=[Edev(1)     0.5*Edev(4) 0.5*Edev(6);
       0.5*Edev(4) Edev(2)     0.5*Edev(5);
       0.5*Edev(6) 0.5*Edev(5) Edev(3)    ]; 
   aux1=0.0;
   aux2=0.0;
for i=1:3
    for j=1:3
        aux1=aux1+Sdev3(i,j)*Sdev3(i,j);
        aux2=aux2+Edev3(i,j)*Edev3(i,j);
    end
end
Ueff=0.5*sqrt(aux1/aux2);
end
end