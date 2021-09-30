function [D1,D2,D3]=colomnedefenhenced(DEP,XI)
%
w1=DEP(1);
w2=DEP(2);
w3=DEP(3);

aux1=XI(3)*(XI(2)-2*XI(1)*XI(2)-XI(2)*XI(2));
aux2=XI(3)*(XI(1)-2*XI(1)*XI(2)-XI(1)*XI(1));
aux3=-2*XI(3);
%            
D1=[w1 w2 0]*aux1;
D2=[w1 w2 0]*aux2;
D3=[0  0 w3]*aux3;
end