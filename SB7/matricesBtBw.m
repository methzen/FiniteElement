%%
%%appel des deux proc�dures pr�c�dentes Bt_Bw
function [Bt,Bw]=matricesBtBw(ni,nj,nk,x1,y1,x2,y2,x3,y3,aire)
      
      x21=x2-x1;
      y21=y2-y1;
      L21=sqrt(x21*x21+y21*y21);
      c4=x21/L21;
      s4=y21/L21;
      x32=x3-x2;
      y32=y3-y2;
      L32=sqrt(x32*x32+y32*y32);
      c5=x32/L32;
      s5=y32/L32;
      x13=x1-x3;
      y13=y1-y3;
      L13=sqrt(x13*x13+y13*y13);
      c6=x13/L13;
      s6=y13/L13;

%%
%%matrice gradient faisant intervenir les rotations au milieu des cot�s 
%%Bt du l'�l�ment DKT6
[Bt]=matriceBt(ni,nj,nk,x21,y21,x32,y32,x13,y13,...
     aire,c4,s4,c5,s5,c6,s6);
      
%% 
%%matrice gradient faisant intervenir les translations normales aux sommets 
%%de l'�l�ment DKT6
[Bw]=matriceBw(aire,c4,s4,c5,s5,c6,s6);

end
