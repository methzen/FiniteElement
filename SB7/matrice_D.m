function [D,C]= matrice_D(E,v)
lambda=v*E/((1-2*v)*(1+v));
mu=E/(2*(1+v));

D=[
    (lambda+2*mu)   lambda          lambda         0       0       0
    lambda          (lambda+2*mu)   lambda         0       0       0
    lambda          lambda          (lambda+2*mu)  0       0       0     
    0               0               0              mu      0       0
    0               0               0              0       mu      0
    0               0               0              0       0       mu  
];

lambdap=E/(1-v*v);
Dp=[
    lambdap   lambdap     0        0       0       0
    lambdap   lambdap     0        0       0       0
    0        0            E        0       0       0     
    0        0            0        mu      0       0
    0        0            0        0       mu      0
    0        0            0        0       0       mu  
];

  C=[ 4/3 -2/3  -2/3  0 0 0;
     -2/3  4/3  -2/3  0 0 0;
     -2/3 -2/3   4/3  0 0 0
     0     0    0   1 0 0;
     0     0    0   0 1 0;
     0     0    0   0 0 1]*mu;
end