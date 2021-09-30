function Benh=Bctenh(XYZ,XI)

[J1,J2,J3]=colomnejacobien(XYZ,XI);
%
aux1=27*(XI(2)-2*XI(1)*XI(2)-XI(2)*XI(2))*XI(3);
aux2=27*(XI(1)-2*XI(1)*XI(2)-XI(1)*XI(1))*XI(3);

%aux1=(2*(XI(1)+XI(2))-1)/2;
%aux2=(2*(XI(1)+XI(2))-1)/2;

aux3=-2*XI(3);

%
Benh=[J1'*aux1;
      J2'*aux2;
      J3'*aux3;
      J1'*aux2+J2'*aux1;
      J1'*aux3+J3'*aux1;
      J2'*aux3+J3'*aux2]; 

% Benh=[aux1 0    0;
%       0    aux2 0;
%       0    0    aux3;
%       aux2 aux1 0 ;
%       0 0 0;
%       0 0 0];  
  
end