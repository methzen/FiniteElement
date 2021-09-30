%%
%%B membrane
%%en input ce sont les translations dans le rep�re global et renvoit l'accroissement 
%%de la d�formation membranaire 
%%n'est pas utilis� pour la stabilisation
function [Bm]=MatriceBm0(xX,xY,xZ,yX,yY,yZ,bx,by)

   Bm(1,1 )=bx(1)*xX; Bm(2,1 )=by(1)*yX; Bm(3,1 )=bx(1)*yX+by(1)*xX;
   Bm(1,2 )=bx(1)*xY; Bm(2,2 )=by(1)*yY; Bm(3,2 )=bx(1)*yY+by(1)*xY;
   Bm(1,3 )=bx(1)*xZ; Bm(2,3 )=by(1)*yZ; Bm(3,3 )=bx(1)*yZ+by(1)*xZ;
   Bm(1,4 )=bx(2)*xX; Bm(2,4 )=by(2)*yX; Bm(3,4 )=bx(2)*yX+by(2)*xX;
   Bm(1,5 )=bx(2)*xY; Bm(2,5 )=by(2)*yY; Bm(3,5 )=bx(2)*yY+by(2)*xY;
   Bm(1,6 )=bx(2)*xZ; Bm(2,6 )=by(2)*yZ; Bm(3,6 )=bx(2)*yZ+by(2)*xZ;
   Bm(1,7 )=bx(3)*xX; Bm(2,7 )=by(3)*yX; Bm(3,7 )=bx(3)*yX+by(3)*xX;
   Bm(1,8 )=bx(3)*xY; Bm(2,8 )=by(3)*yY; Bm(3,8 )=bx(3)*yY+by(3)*xY;
   Bm(1,9 )=bx(3)*xZ; Bm(2,9 )=by(3)*yZ; Bm(3,9 )=bx(3)*yZ+by(3)*xZ;
   Bm(1,10)=bx(4)*xX; Bm(2,10)=by(4)*yX; Bm(3,10)=bx(4)*yX+by(4)*xX;
   Bm(1,11)=bx(4)*xY; Bm(2,11)=by(4)*yY; Bm(3,11)=bx(4)*yY+by(4)*xY;
   Bm(1,12)=bx(4)*xZ; Bm(2,12)=by(4)*yZ; Bm(3,12)=bx(4)*yZ+by(4)*xZ;
   Bm(1,13)=bx(5)*xX; Bm(2,13)=by(5)*yX; Bm(3,13)=bx(5)*yX+by(5)*xX;
   Bm(1,14)=bx(5)*xY; Bm(2,14)=by(5)*yY; Bm(3,14)=bx(5)*yY+by(5)*xY;
   Bm(1,15)=bx(5)*xZ; Bm(2,15)=by(5)*yZ; Bm(3,15)=bx(5)*yZ+by(5)*xZ;
   Bm(1,16)=bx(6)*xX; Bm(2,16)=by(6)*yX; Bm(3,16)=bx(6)*yX+by(6)*xX;
   Bm(1,17)=bx(6)*xY; Bm(2,17)=by(6)*yY; Bm(3,17)=bx(6)*yY+by(6)*xY;
   Bm(1,18)=bx(6)*xZ; Bm(2,18)=by(6)*yZ; Bm(3,18)=bx(6)*yZ+by(6)*xZ;

end