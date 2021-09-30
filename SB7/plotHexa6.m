function plotHexa6(XYZ,c)
X=XYZ(:,1);
Y=XYZ(:,2);
Z=XYZ(:,3);

plot3([X(1),X(2)],[Y(1),Y(2)],[Z(1),Z(2)],c)
hold on % pour ne pas effacer ce qui a �t� affich� pr�c�demment
plot3([X(2),X(3)],[Y(2),Y(3)],[Z(2),Z(3)],c)
plot3([X(3),X(1)],[Y(3),Y(1)],[Z(3),Z(1)],c)
plot3([X(4),X(5)],[Y(4),Y(5)],[Z(4),Z(5)],c)
plot3([X(5),X(6)],[Y(5),Y(6)],[Z(5),Z(6)],c)
plot3([X(6),X(4)],[Y(6),Y(4)],[Z(6),Z(4)],c)
plot3([X(1),X(4)],[Y(1),Y(4)],[Z(1),Z(4)],c)
plot3([X(2),X(5)],[Y(2),Y(5)],[Z(2),Z(5)],c)
plot3([X(3),X(6)],[Y(3),Y(6)],[Z(3),Z(6)],c)


% pour visualiser avec les bonnes dimensions
 view(3)
 set(gca,'DataAspectRatio',[1 1 1])
 set(gca,'visible','on')
 set(gcf,'Color','w')

end
