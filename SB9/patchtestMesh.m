
C=[1 2 3 4
    5 6 2 1
    6 7 3 2
    7 8 4 3
    8 5 1 4];
C1=C+1;

N=[0.04 0.02 0.0
    0.18 0.03 0.0
    0.16 0.08 0.0
    0.08 0.08 0.0
    0.0  0.0  0.0
    0.24 0.0  0.0
    0.24 0.12 0.0
    0.0  0.12 0.0];
N1=N;
N1(:,3)=0.0005;

LE=[C C1];
XYZ=[N;N1];