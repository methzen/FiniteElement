function BB=calculmatBSB6(XYZ,XI)
r=XI(1);
s=XI(2);
zeta=XI(3);
Jac=Jacobien(XYZ,XI);
Nr=[1-zeta,0,-(1-zeta),1+zeta,0,-(1+zeta)]/2;
Ns=[0,1-zeta,-(1-zeta),0, 1+zeta,-(1+zeta)]/2;
Nzeta=[-r,-s,-(1-r-s),r,s,(1-r-s)]/2;

N=inv(Jac)*[Nr;Ns;Nzeta];
Nx=N(1,:);
Ny=N(2,:);
Nz=N(3,:);
zero6=zeros(1,6);
BB=zeros(6,18);
for i=1:6
col=3*(i-1)+1:3*(i-1)+3;
BB(:,col)=[Nx(i) zero6(i) zero6(i);
           zero6(i) Ny(i) zero6(i);
           zero6(i) zero6(i) Nz(i);
           Ny(i) Nx(i) zero6(i);
           zero6(i) Nz(i) Ny(i);           
           Nz(i) zero6(i) Nx(i)];
end
    
