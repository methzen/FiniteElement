function T=matrice_trans_rot(h1,h2,h3,...
            V1xx,V1xy,V1xz,V1yx,V1yy,V1yz,V1zx,V1zy,V1zz,...
            V2xx,V2xy,V2xz,V2yx,V2yy,V2yz,V2zx,V2zy,V2zz,...
            V3xx,V3xy,V3xz,V3yx,V3yy,V3yz,V3zx,V3zy,V3zz)
T=zeros(15,18);
%%       
T(1,1)=1.0;
T(4,1)=V1yx/h1;
T(5,1)=-V1xx/h1;

T(4,2)=V1yy/h1;
T(5,2)=-V1xy/h1;
T(6,2)=1.0;

T(4,3)=V1yz/h1;
T(5,3)=-V1xz/h1;
T(11,3)=1.0;
%%
T(2,4)=1.0;
T(9,4)=V2yx/h2;
T(10,4)=-V2xx/h2;

T(7,5)=1;
T(9,5)=V2yy/h2;
T(10,5)=-V2xy/h2;

T(9,6)=V2yz/h2;
T(10,6)=-V2xz/h2;
T(12,6)=1.0;
%%
%
T(3,7)=1.0;
T(14,7)=V3yx/h3;
T(15,7)=-V3xx/h3;

T(8,8)=1;
T(14,8)=V3yy/h3;
T(15,8)=-V3xy/h3;

T(13,9)=1.0;
T(14,9)=V3yz/h3;
T(15,9)=-V3xz/h3;
%%
T(1,1+9)=1.0;
T(4,1+9)=-V1yx/h1;
T(5,1+9)=V1xx/h1;

T(4,2+9)=-V1yy/h1;
T(5,2+9)=V1xx/h1;
T(6,2+9)=1.0;

T(4,3+9)=-V1yz/h1;
T(5,3+9)=V1xz/h1;
T(11,3+9)=1.0;
%%
T(2,4+9)=1.0;
T(9,4+9)=-V2yx/h2;
T(10,4+9)=V2xx/h2;

T(7,5+9)=1;
T(9,5+9)=-V2yy/h2;
T(10,5+9)=V2xy/h2;

T(9,6+9)=-V2yz/h2;
T(10,6+9)=V2xz/h2;
T(12,6+9)=1.0;
%%
T(3,7+9)=1.0;
T(14,7+9)=-V3yx/h3;
T(15,7+9)=V3xx/h3;

T(8,8+9)=1;
T(14,8+9)=-V3yy/h3;
T(15,8+9)=V3xy/h3;

T(13,9+9)=1.0;
T(14,9+9)=-V3yz/h3;
T(15,9+9)=V3xz/h3;
%%
 
T=0.5*T;
end