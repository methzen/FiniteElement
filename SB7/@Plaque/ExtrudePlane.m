function [ Nodes3D,Connect3D] = ExtrudePlane( obj)
%Extrusion est une fonction qui renvoie une table de connectivit� d'un
%maillage 3D � partir d'un maillage plan.
%   Detailed explanation goes here
Nodes2D=obj.Nodes2D;
Con=obj.Connect2D;
e=obj.Epaiss;
NbrEtage=obj.discretz;

h=e/2;
p=size(Con,2);
M1=Con(:,2:p);
inc=max(max(M1));
M2=M1+inc;


[NbrNoeud2D,n]=size(Nodes2D);

NbrNoeudTotal=NbrNoeud2D*(NbrEtage+1);
Nodes3D=sparse(NbrNoeudTotal,n);
Nodes3D(:,1)=1:NbrNoeudTotal;
Nodes2D(:,4)=-h;
for j=1:(NbrEtage+1)
    
   a=NbrNoeud2D*(j-1)+1; b=j*NbrNoeud2D;
   Nodes3D(a:b,2:3)=Nodes2D(:,2:3);
   Nodes3D(a:b,4)=Nodes2D(:,4)+e*(j-1)/NbrEtage;
end    


E=[M1 M2];

[p,q]=size(E);
Connect3D=zeros(NbrEtage*p,q);

for i=1:NbrEtage
    a=p*(i-1)+1; b=i*p;
    Connect3D(a:b,:)=E+(i-1)*inc;

end

