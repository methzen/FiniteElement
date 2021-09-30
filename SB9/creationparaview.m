clear all
clc
% [XYZ,LE,Cx,Cy]=hemisphere;
% DISPTD=deplhemis;
%  NOPRVW=0;
%  U2=zeros(10,1);
%  U50=zeros(10,1); 
%  for i=1:10
% A=DISPTD(1:578,2:4,i);
% U2(i)=A(2,1);
% U50(i)=A(50,2);
% % B=A';
% %  NOPRVW=paraview_output(LE,XYZ,B,NOPRVW);
%  end

% A=zeros(100,1);
% B=zeros(100,1);
% C=zeros(100,1);
%  U=deplcylindre;
%  for i=1:100
%      x=U(58,2,i);
%      y=U(58,3,i);
% A(i)=sqrt(x*x+y*y+z*z);
%      x=U(2,2,i);
%      y=U(2,3,i);
%      z=U(2,4,i);
% B(i)=sqrt(x*x+y*y+z*z);
%      x=U(3,2,i);
%      y=U(3,3,i);
%      z=U(3,4,i);
% C(i)=sqrt(x*x+y*y+z*z);
%  end

A=zeros(100,1);
B=zeros(100,1);
C=zeros(100,1);
 U=deplcylindre;
 for i=1:100
     y=U(58,3,i);
A(i)=abs(y);
     x=U(2,2,i);
B(i)=sqrt(x*x);
     x=U(3,2,i);
C(i)=sqrt(x*x);
 end
D=0:0.4:40;
M=[0 0 0;A B C]