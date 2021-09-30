function [Nodes2D,Connect2D,Numerot2D] = DiscretPlane(obj)
%UNTITLED2 Summary of this function goes here
% Lx Dimension de la plaque selon x
% Ly diemension de la plaque selon y
% x est le nombre d'élément selon x
% y est le nombre d'élément selon y
% if isinteger(x)==1
Lx=obj.Long;
Ly=obj.Larg;
x=obj.discretx;
y=obj.discrety;

NbrElem=x*y;
NbrNoeu=(x+1)*(y+1);    
M=1:NbrNoeu ;

A=reshape(M,x+1,y+1);
Numerot2D=A';

Nodes2D=sparse(NbrNoeu,4); %tableau des coordonnées des noeuds

Con=zeros(NbrElem,5);
boucle=1;
for j=1:(y+1)
     for i=1:(x+1)
        m=Numerot2D(j,i);
        Y=-Ly/2+(j-1)*Ly/y;
        X=-Lx/2+(i-1)*Lx/x;
        Nodes2D(m,:)=[m X Y 0];
    end
end

while boucle <= NbrElem % boucle creant la table de connectivité

for j=1:y
        for i=1:x
            N1=Numerot2D(j,i);
            N2=Numerot2D(j,i)+1;
            N3=N2+(x+1);
            N4=N1+(x+1);
            Con(boucle,:)=[boucle N1 N2 N3 N4];
            boucle=boucle+1;
        end
end
    
    
end

Connect2D=Con;