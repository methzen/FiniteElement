    function [E13,E23]=defct(XYZ,DSP,XI)
              [J1,J2,J3]=colomnejacobien(XYZ,XI);
              [D1,D2,D3]=colomneDefconv(DSP,XI); 
              E13=dot(J1,D3)+dot(J3,D1)+dot(D1,D3);
              E23=dot(J2,D3)+dot(J3,D2)+dot(D2,D3);
    end