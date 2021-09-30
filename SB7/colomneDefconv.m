function [D1,D2,D3]=colomneDefconv(DSP,XI)
%
[D10,D1ZETA,...
 D20,D2ZETA,...
 D30,D3ETA,D3XI]=deformaconva(DSP);
%           
D1=D10+XI(3)*D1ZETA;
D2=D20+XI(3)*D2ZETA;
D3=D30+XI(2)*D3ETA+XI(1)*D3XI;
%
D1=D1';
D2=D2';
D3=D3';
end