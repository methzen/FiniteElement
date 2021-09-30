function [Cprisme,XYZ]=Mailcubeenprisme(Ccube,XYZ,aux)

if aux==1
Cprisme=zeros(2*size(Ccube,1),6);
%
for i=1:size(Ccube,1)
    k=2*i-1:2*i;
    f=Ccube(i,:);
    Cprisme(k,:)=[f(3) f(1) f(2) f(7) f(5) f(6);
                  f(1) f(3) f(4) f(5) f(7) f(8)];
end
%
elseif aux==4 
ma=max(max(Ccube));
nomb=2*size(Ccube,1);
Noesup=ma+1:nomb+ma;
Cprisme=zeros(4*size(Ccube,1),6);
for i=1:size(Ccube,1)
    c=2*(i-1)+1;
    nnin=Noesup(c);
    nnsu=Noesup(c+1);
    k=4*(i-1)+1:4*i;
    f=Ccube(i,:);
    Cprisme(k,:)=[f(1) nnin f(4) f(5) nnsu f(8);
                  f(2) nnin f(1) f(6) nnsu f(5);
                  f(3) nnin f(2) f(7) nnsu f(6);
                  f(4) nnin f(3) f(8) nnsu f(7)];
end
[m,n]=size(XYZ);
XYZ1=zeros(m+2*size(Ccube,1),n);
XYZ1(1:m,1:n)=XYZ;
for i=1:size(Ccube,1)
    f=Ccube(i,:);
        c=2*(i-1)+1;
        nnin=Noesup(c);
        nnsu=Noesup(c+1);
        XYZ1(nnin,1)=sum(XYZ(f(1:4),1))/4;
        XYZ1(nnin,2)=sum(XYZ(f(1:4),2))/4;
        XYZ1(nnin,3)=sum(XYZ(f(1:4),3))/4;
        XYZ1(nnsu,1)=sum(XYZ(f(5:8),1))/4;
        XYZ1(nnsu,2)=sum(XYZ(f(5:8),2))/4;
        XYZ1(nnsu,3)=sum(XYZ(f(5:8),3))/4;        
end
XYZ=XYZ1;
end
end