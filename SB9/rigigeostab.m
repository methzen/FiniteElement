function Kgstab=rigigeostab(SH1,SH2,SH23,SH13,BGxi,BGeta,BGxizeta,BGetazeta)

%Coordonnées initiale en entrée
%%
       SIG1=[SH1(1) SH1(4) SH1(6);
             SH1(4) SH1(2) SH1(5);
             SH1(6) SH1(5) SH1(3)];
        SHEAD1=zeros(9);
        SHEAD1(1:3,1:3)=SIG1;
        SHEAD1(4:6,4:6)=SIG1;
        SHEAD1(7:9,7:9)=SIG1;
%%        
       SIG2=[SH2(1) SH2(4) SH2(6);
             SH2(4) SH2(2) SH2(5);
             SH2(6) SH2(5) SH2(3)];
        SHEAD2=zeros(9);
        SHEAD2(1:3,1:3)=SIG2;
        SHEAD2(4:6,4:6)=SIG2;
        SHEAD2(7:9,7:9)=SIG2;
%%        
  SIG23=[SH23(1) SH23(4) SH23(6);
         SH23(4) SH23(2) SH23(5);
         SH23(6) SH23(5) SH23(3)];
         SHEAD23=zeros(9);
         SHEAD23(1:3,1:3)=SIG23;
         SHEAD23(4:6,4:6)=SIG23;
         SHEAD23(7:9,7:9)=SIG23;
%%        
   SIG13=[SH13(1) SH13(4) SH13(6);
         SH13(4) SH13(2) SH13(5);
         SH13(6) SH13(5) SH13(3)];
         SHEAD13=zeros(9);
         SHEAD13(1:3,1:3)=SIG13;
         SHEAD13(4:6,4:6)=SIG13;
         SHEAD13(7:9,7:9)=SIG13;
%%
 Kgstab=8*(BGxi'*SHEAD1*BGxi+BGeta'*SHEAD2*BGeta)/3+...
        8*(BGetazeta'*SHEAD23*BGetazeta+BGxizeta'*SHEAD13*BGxizeta)/9.0;        
end
