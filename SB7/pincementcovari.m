    function E33=pincementcovari(XYZ,DSP,XI)
             [~,~,J3]=colomnejacobien(XYZ,XI);
             [~,~,D3]=colomneDefconv(DSP,XI); 
             E33=dot(J3,D3)+0.5*dot(D3,D3);
    end