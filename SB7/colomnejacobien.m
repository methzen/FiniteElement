function [J1,J2,J3]=colomnejacobien(XYZ,XI)
%


% Jac=Jacobien(XYZ,XI);
% 
% J1=Jac(:,1);
% J2=Jac(:,2);
% J3=Jac(:,3);

[J10,J1ZETA,...
 J20,J2ZETA,...
 J30,J3ETA,J3XI]=jacobienDecompose(XYZ);
%            
J1=J10+XI(3)*J1ZETA;
J2=J20+XI(3)*J2ZETA;
J3=J30+XI(2)*J3ETA+XI(1)*J3XI;
%
J1=J1';
J2=J2';
J3=J3';
end