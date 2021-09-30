%%
%%parail pour le pincement
function [Bpc]=MatriceBp0(zX,zY,zZ,bz)
      Bpc(1 )=bz(1)*zX;
      Bpc(2 )=bz(1)*zY;
      Bpc(3 )=bz(1)*zZ;
      Bpc(4 )=bz(2)*zX;
      Bpc(5 )=bz(2)*zY;
      Bpc(6 )=bz(2)*zZ;
      Bpc(7 )=bz(3)*zX;
      Bpc(8 )=bz(3)*zY;
      Bpc(9 )=bz(3)*zZ;
      Bpc(10)=bz(4)*zX;
      Bpc(11)=bz(4)*zY;
      Bpc(12)=bz(4)*zZ;
      Bpc(13)=bz(5)*zX;
      Bpc(14)=bz(5)*zY;
      Bpc(15)=bz(5)*zZ;
      Bpc(16)=bz(6)*zX;
      Bpc(17)=bz(6)*zY;
      Bpc(18)=bz(6)*zZ;
end