function DFF=dshapePrisme(XI)
r=XI(1);
s=XI(2);
t=XI(3);
% DFF=[(1-t) 0 -(1-t) (1+t) 0  -(1+t);
%       0 (1-t) -(1-t) 0 (1+t) -(1+t);
%      -r  -s -(1-r-s) r s (1-r-s)]/2;

 
DFF=[1-t,0,-(1-t),1+t,0,-(1+t);
     0,1-t,-(1-t),0, 1+t,-(1+t);
    -r,-s,-(1-r-s),r,s,(1-r-s)]/2;
end