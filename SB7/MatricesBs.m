%%
%%matrice gradient pour la stabilisation
%%Bsx pour q1x,q1y; Bsy pour q1y,q2y; Bsz pour q1z,q2z "d�formations" de stabilisation
function [Bsx,Bsy,Bsz]=MatricesBs(xX,xY,xZ,yX,yY,yZ,zX,zY,zZ,Vgamma)
  for j=1:6
 %%Bsx pour q1x,q1y; Bsy pour q1y,q2y; Bsz pour q1z,q2z
       k=3*(j-1);
       Bsx(1,k+1)=xX*Vgamma(j,1); Bsx(1,k+2)=xY*Vgamma(j,1); Bsx(1,k+3)=xZ*Vgamma(j,1);
       Bsx(2,k+1)=xX*Vgamma(j,2); Bsx(2,k+2)=xY*Vgamma(j,2); Bsx(2,k+3)=xZ*Vgamma(j,2);

       Bsy(1,k+1)=yX*Vgamma(j,1); Bsy(1,k+2)=yY*Vgamma(j,1); Bsy(1,k+3)=yZ*Vgamma(j,1);
       Bsy(2,k+1)=yX*Vgamma(j,2); Bsy(2,k+2)=yY*Vgamma(j,2); Bsy(2,k+3)=yZ*Vgamma(j,2);

       Bsz(1,k+1)=zX*Vgamma(j,1); Bsz(1,k+2)=zY*Vgamma(j,1); Bsz(1,k+3)=zZ*Vgamma(j,1);
       Bsz(2,k+1)=zX*Vgamma(j,2); Bsz(2,k+2)=zY*Vgamma(j,2); Bsz(2,k+3)=zZ*Vgamma(j,2);
  end
end