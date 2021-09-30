function plotHexa8(XYZ,c)
X=XYZ(:,1);
Y=XYZ(:,2);
Z=XYZ(:,3);

plot3([X(1),X(2)],[Y(1),Y(2)],[Z(1),Z(2)],c)
hold on % pour ne pas effacer ce qui a été affiché précédemment
plot3([X(2),X(3)],[Y(2),Y(3)],[Z(2),Z(3)],c)
plot3([X(3),X(4)],[Y(3),Y(4)],[Z(3),Z(4)],c)
plot3([X(4),X(1)],[Y(4),Y(1)],[Z(4),Z(1)],c)
plot3([X(5),X(6)],[Y(5),Y(6)],[Z(5),Z(6)],c)
plot3([X(6),X(7)],[Y(6),Y(7)],[Z(6),Z(7)],c)
plot3([X(7),X(8)],[Y(7),Y(8)],[Z(7),Z(8)],c)
plot3([X(8),X(5)],[Y(8),Y(5)],[Z(8),Z(5)],c)
plot3([X(1),X(5)],[Y(1),Y(5)],[Z(1),Z(5)],c)
plot3([X(2),X(6)],[Y(2),Y(6)],[Z(2),Z(6)],c)
plot3([X(3),X(7)],[Y(3),Y(7)],[Z(3),Z(7)],c)
plot3([X(4),X(8)],[Y(4),Y(8)],[Z(4),Z(8)],c)


% pour visualiser avec les bonnes dimensions
 view(3)
 set(gca,'DataAspectRatio',[1 1 1])
set(gca,'visible','on')
 set(gcf,'Color','w')

end
