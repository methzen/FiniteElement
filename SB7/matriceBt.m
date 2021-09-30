%%
%%matrice gradient faisant intervenir les rotations au milieu des cot�s 
%%Bt du l'�l�ment DKT6
function [Bt]=matriceBt(ni,nj,nk,x21,y21,x32,y32,x13,y13,...
          aire,c4,s4,c5,s5,c6,s6)
      
      Bt(1,1)=(y21)*s4/aire;
      Bt(1,2)=(y32)*s5/aire;
      Bt(1,3)=(y13)*s6/aire;
      Bt(2,1)=(x21)*c4/aire;
      Bt(2,2)=(x32)*c5/aire;
      Bt(2,3)=(x13)*c6/aire;
      Bt(3,1)=((-y21)*c4+(-x21)*s4)/aire;
      Bt(3,2)=((-y32)*c5+(-x32)*s5)/aire;
      Bt(3,3)=((-y13)*c6+(-x13)*s6)/aire;
      if(ni>nj)
        Bt(1,1)=-Bt(1,1);
        Bt(2,1)=-Bt(2,1);
        Bt(3,1)=-Bt(3,1);
      end
      if(nj>nk)
        Bt(1,2)=-Bt(1,2);
        Bt(2,2)=-Bt(2,2);
        Bt(3,2)=-Bt(3,2);
      end
      if(nk>ni)
        Bt(1,3)=-Bt(1,3);
        Bt(2,3)=-Bt(2,3);
        Bt(3,3)=-Bt(3,3);
      end
end