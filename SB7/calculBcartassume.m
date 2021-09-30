function B=calculBcartassume(XYZ,xyz,XI)
%
B=calculB(xyz,XI);
%% cisaillement
Bct=Bchassumect(xyz,XI);
%Bct=Bctboisse1(xyz,XI);
%Bct=MatriceBc0(XYZ)*(1-XI(3)*XI(3))*(5/6);
B(5:6,:)=Bct;
%% pinching
Bp=Bpinch(XYZ,XI);
B(3,:)=Bp;
%% membrane
%Bm=Bmembrane(XYZ,XI);
%B(1:2,:)=Bm(1:2,:);
%B(4,:)=Bm(3,:);
%%
%T=matriceT(XYZ,XI);
T0=matriceT(XYZ,[0 0 0]);
%
T0(5:6,1:4)=0;
T0(1:4,5:6)=0;
B7=T0*[0 0 XI(3) 0 0 0]';
B=[T0*B B7];
end