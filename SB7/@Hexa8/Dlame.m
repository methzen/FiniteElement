function D= Dlame(obj,E,v)

lambda=v*E/((1-2*v)*(1+v));
mu=E/(2*(1+v));
D=[
    (lambda+2*mu)   lambda          lambda          0       0       0
    lambda          (lambda+2*mu)   lambda          0       0       0
    lambda          lambda          (lambda+2*mu)   0       0       0
    0               0               0               mu      0       0
    0               0               0               0       mu      0
    0               0               0               0       0       mu  
];
end