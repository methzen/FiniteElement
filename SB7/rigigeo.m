function Kg=rigigeo(S,XYZ,XI)


% Cette matrice Kg est celle sur un points de Gauss
T=matriceT(XYZ,XI);
DFF=dshapePrisme(XI);
%%
Kg=zeros(18,18);
for I=1:6
    for J=1:6
    LIN=3*(I-1)+1:3*(I-1)+3;
    COL=3*(J-1)+1:3*(J-1)+3;
    B=[DFF(1,I)*DFF(1,J);
       DFF(2,I)*DFF(2,J);
       DFF(3,I)*DFF(3,J);
       DFF(1,I)*DFF(2,J)+DFF(2,I)*DFF(1,J);
       DFF(3,I)*DFF(2,J)+DFF(2,I)*DFF(3,J);
       DFF(1,I)*DFF(3,J)+DFF(3,I)*DFF(1,J)];
    Bg=T*B;
       Kg(LIN,COL)=dot(Bg,S)*eye(3);
    end
end

end
