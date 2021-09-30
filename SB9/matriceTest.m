function T=matriceTest(J)

J11=J(1,1);
J12=J(1,2);
J13=J(1,3);
J21=J(2,1);
J22=J(2,2);
J23=J(2,3);
J31=J(3,1);
J32=J(3,2);
J33=J(3,3);

T=[  J11*J11   J21*J21   J31*J31         J11*J21         J21*J31         J11*J31         
     J12*J12   J22*J22   J32*J32         J12*J22         J22*J32         J12*J32         
     J13*J13   J23*J23   J33*J33         J13*J23         J23*J33         J13*J33         
   2*J11*J12 2*J21*J22 2*J31*J32 J11*J22+J21*J12 J22*J31+J21*J32 J11*J32+J31*J12
   2*J12*J13 2*J22*J23 2*J32*J33 J13*J22+J12*J23 J22*J33+J32*J23 J12*J33+J32*J13
   2*J11*J13 2*J21*J23 2*J31*J33 J13*J21+J11*J23 J23*J31+J21*J33 J13*J31+J11*J33];
end