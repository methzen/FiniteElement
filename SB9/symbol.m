gI1=sym('gI1');
gI2=sym('gI2');
gI3=sym('gI3');
hI1=sym('hI1');
hI2=sym('hI2');
hI3=sym('hI3');
hI4=sym('hI4');

gJ1=sym('gJ1');
gJ2=sym('gJ2');
gJ3=sym('gJ3');
hJ1=sym('hJ1');
hJ2=sym('hJ2');
hJ3=sym('hJ3');
hJ4=sym('hJ4');
zeta=sym('zeta');
eta=sym('eta');
xi=sym('xi');

NI1=gI1+eta*hI1+zeta*hI3+eta*zeta*hI4;
NI2=gI2+xi*hI1+zeta*hI2+xi*zeta*hI4;
NI3=gI3+eta*hI2+xi*hI3+eta*xi*hI4;

NJ1=gJ1+eta*hJ1+zeta*hJ3+eta*zeta*hJ4;
NJ2=gJ2+xi*hJ1+zeta*hJ2+xi*zeta*hJ4;
NJ3=gJ3+eta*hJ2+xi*hJ3+eta*xi*hJ4;

G11=NI1*NJ1
G22=NI2*NJ2
G33=NI3*NJ3
G12=NI2*NJ1+NI1*NJ2
G23=NI3*NJ2+NI2*NJ3
G13=NI3*NJ1+NI1*NJ3

