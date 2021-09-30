function [Ec0,Eczeta,Eczz,EcXI,EcETA,EcXIETA,EcETAZETA,EcXIZETA]=deformatEc(DSP,XYZ)

% iCI LES COORDONNÉES SONT LES COORDONNÉES INITIALES
%  ALERTE ECRITURE VERIFIEE
%
[J10,J1ETA,J1ZETA,J1ETAZETA,...
 J20,J2XI,J2ZETA,J2XIZETA,...
 J30,J3ETA,J3XI,J3XIETA]=jacobienDecompose(XYZ);

 [D10,D1ETA,D1ZETA,D1ETAZETA,...
  D20,D2XI,D2ZETA,D2XIZETA,...
  D30,D3ETA,D3XI,D3XIETA]=deformaconva(DSP);
%%
% A à D
XIAD=[-1 -1 0
       1 -1 0
       1  1 0
      -1  1 0];
% E TO H
XIEH=[-1  0 -1 
       1  0 -1
       1  0  1
      -1  0  1];
% J TO M
XIJM=[0 -1 -1
      0  1 -1
      0  1  1
      0 -1  1];
%%
Ec0=zeros(6,1);
Ec0(1)=dot(J10,D10)+0.5*dot(D10,D10);
Ec0(2)=dot(J20,D20)+0.5*dot(D20,D20);
Ec0(4)=dot(J10,D20)+dot(J20,D10)+dot(D10,D20);
%%
Eczeta=zeros(6,1);
Eczeta(1)=dot(J10,D1ZETA)+dot(D10,J1ZETA)+dot(D10,D1ZETA);
Eczeta(2)=dot(J20,D2ZETA)+dot(D20,J2ZETA)+dot(D20,D2ZETA);
Eczeta(4)=dot(J10,D2ZETA)+dot(J1ZETA,D20)+dot(J20,D1ZETA)+dot(J2ZETA,D10)...
          +dot(D10,D2ZETA)+dot(D1ZETA,D20);

%%
Eczz=zeros(6,1);
Eczz(1)=dot(J1ZETA,D1ZETA)+0.5*dot(D1ZETA,D1ZETA);
Eczz(2)=dot(J2ZETA,D2ZETA)+0.5*dot(D2ZETA,D2ZETA);
Eczz(4)=dot(J1ZETA,D2ZETA)+dot(J2ZETA,D1ZETA)+dot(D2ZETA,D1ZETA);
%% STABILISATION XI
EcXI=zeros(6,1);
EcXI(2)=dot(J20,D2XI)+dot(J2XI,D20)+dot(D2XI,D20);
EcXI(4)=dot(J10,D2XI)+dot(J2XI,D10)+dot(D2XI,D10);
%% STABILISATION ETA
EcETA=zeros(6,1);
EcETA(1)=dot(J10,D1ETA)+dot(J1ETA,D10)+dot(D1ETA,D10);
EcETA(4)=dot(J20,D1ETA)+dot(J1ETA,D20)+dot(D1ETA,D20);
%% STABILISATION XIETA
EcXIETA=zeros(6,1);
EcXIETA(4)=dot(J1ETA,D2XI)+dot(J2XI,D1ETA)+dot(D2XI,D1ETA);  
%% STABILISATION ETAZETA
EcETAZETA=zeros(6,1);
EcETAZETA(1)=dot(J10,D1ETAZETA)+dot(J1ETA,D1ZETA)+dot(J1ZETA,D1ETA)...
         +dot(J1ETAZETA,D10)+dot(D1ETAZETA,D10)+dot(D1ETA,D1ZETA);
EcETAZETA(4)=dot(J20,D1ETAZETA)+dot(J1ETA,D2ZETA)+dot(J2ZETA,D1ETA)...
         +dot(J1ETAZETA,D20)+dot(D1ETAZETA,D20)+dot(D1ETA,D2ZETA);
%% STABILISATION XIZETA
EcXIZETA=zeros(6,1);
EcXIZETA(2)=dot(J20,D2XIZETA)+dot(J2XI,D2ZETA)+dot(J2ZETA,D2XI)...
         +dot(J2XIZETA,D20)+dot(D2XIZETA,D20)+dot(D2XI,D2ZETA);
EcXIZETA(4)=dot(J10,D2XIZETA)+dot(J2XI,D1ZETA)+dot(J1ZETA,D2XI)...
         +dot(J2XIZETA,D10)+dot(D2XIZETA,D10)+dot(D2XI,D1ZETA);     
%% CHAMPS ASSUMES
for AD=1:4
   [~,~,J3]=colomnejacobien(XYZ,XIAD(AD,:));
   [~,~,D3]=colomneDefconv(DSP,XIAD(AD,:));
   Ec0(3)= Ec0(3)+(dot(J3,D3)+0.5*dot(D3,D3))/4;
   EcXI(3)= EcXI(3)+XIAD(AD,1)*(dot(J3,D3)+0.5*dot(D3,D3))/4;
   EcETA(3)= EcETA(3)+XIAD(AD,2)*(dot(J3,D3)+0.5*dot(D3,D3))/4;
   EcXIETA(3)= EcXIETA(3)+XIAD(AD,2)*XIAD(AD,1)*(dot(J3,D3)+0.5*dot(D3,D3))/4;
end

for EH=1:4
   [~,J2,J3]=colomnejacobien(XYZ,XIEH(EH,:));
   [~,D2,D3]=colomneDefconv(DSP,XIEH(EH,:));
   Ec0(5)= Ec0(5)+(dot(J3,D2)+dot(J2,D3)+dot(D3,D2))/4;
   Eczeta(5)= Eczeta(5)+XIEH(EH,3)*(dot(J3,D2)+dot(J2,D3)+dot(D3,D2))/4;
   EcXI(5)= EcXI(5)+XIEH(EH,1)*(dot(J2,D3)+dot(J3,D2)+dot(D2,D3))/4;
   EcXIZETA(5)= EcXIZETA(5)+XIEH(EH,1)*XIEH(EH,3)*(dot(J2,D3)+dot(J3,D2)+dot(D2,D3))/4;
end

for JM=1:4
   [J1,~,J3]=colomnejacobien(XYZ,XIJM(JM,:));
   [D1,~,D3]=colomneDefconv(DSP,XIJM(JM,:));
   Ec0(6)= Ec0(6)+(dot(J3,D1)+dot(J1,D3)+dot(D3,D1))/4;
   Eczeta(6)= Eczeta(6)+XIJM(JM,3)*(dot(J3,D1)+dot(J1,D3)+dot(D3,D1))/4;
   EcETA(6)= EcETA(6)+XIJM(JM,2)*(dot(J1,D3)+dot(J3,D1)+dot(D1,D3))/4;
   EcETAZETA(6)= EcETAZETA(6)+XIJM(JM,2)*XIJM(JM,3)*(dot(J1,D3)+dot(J3,D1)+dot(D1,D3))/4;
end

end