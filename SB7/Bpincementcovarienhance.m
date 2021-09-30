function B33=Bpincementcovarienhance(XYZ,XI)
     [~,~,J3]=colomnejacobien(XYZ,XI);
     B33=J3*(-2*XI(3));
end