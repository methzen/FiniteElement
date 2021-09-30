function [rho,iter]=linesearchtest()

rho0=.50000001;
rho1=3;
g1=(rho1-0.5)^2;
g0=(rho0-0.5)^2;

iter=0;
while g1>g0
    
    alpha=(rho1-rho0)/(g1-g0);
    rho=rho1-alpha*g1;
    g1=(rho-0.5)^2;
    rho0=rho1;
    rho1=rho;
    iter=iter+1;
end


end