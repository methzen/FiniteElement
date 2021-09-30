function B33I=Bpincementcovari(XYZ,XI,I)
      g3=[0 0 -1 0 0 1]/2;
      h1=[-1 0 1 1 0 -1]/2;
      h2=[0 -1 1 0 1 -1]/2;
      [~,~,J3]=colomnejacobien(XYZ,XI);
      B33I=J3*(g3(I)+h1(I)*XI(1)+h2(I)*XI(2));
end