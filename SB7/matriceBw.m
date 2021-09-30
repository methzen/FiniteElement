%% 
%%matrice gradient faisant intervenir les translations normales aux sommets 
%%de l'�l�ment DKT6
function [Bw]=matriceBw(aire,c4,s4,c5,s5,c6,s6)

      Bw(1,1)=(s4*c4-s6*c6)/aire;
      Bw(2,1)=-Bw(1,1);
      Bw(1,2)=(s5*c5-s4*c4)/aire;
      Bw(2,2)=-Bw(1,2);
      Bw(1,3)=(s6*c6-s5*c5)/aire;
      Bw(2,3)=-Bw(1,3);
      Bw(3,1)=(s4*s4-c4*c4+c6*c6-s6*s6)/aire;
      Bw(3,2)=(s5*s5-c5*c5+c4*c4-s4*s4)/aire;
      Bw(3,3)=(s6*s6-c6*c6+c5*c5-s5*s5)/aire;

end