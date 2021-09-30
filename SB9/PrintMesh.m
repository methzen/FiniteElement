function PrintMesh(Nu,XYZ, Ec,scale)

                [m,n]=size(Ec);
                coord=zeros(8,3);
               % figure('Maillage')
                figure('Name','Maillage','NumberTitle','off')
                set(gcf,'color','white')
                                %% maillage initiale
                 for i=1:m
                    for j=1:n
                        k=Ec(i,j);
                        %u=(k-1)*3+1;
                        coord(j,1)=XYZ(k,1);
                        coord(j,2)=XYZ(k,2);
                        coord(j,3)=XYZ(k,3);
                    end
                    X= coord(:,1)';
                    Y= coord(:,2)';
                    Z= coord(:,3)';

                    plotHexa8(X,Y,Z,'b')
                    set(gca,'Visible','on');
%                     set(gca,'xtick',[],'ytick',[])
                end
                for i=1:m
                    for j=1:n
                        k=Ec(i,j);
                        u=(k-1)*3+1;
                        coord(j,1)=scale*Nu(u)+XYZ(k,1);
                        coord(j,2)=scale*Nu(u+1)+XYZ(k,2);
                        coord(j,3)=scale*Nu(u+2)+XYZ(k,3);
                    end
                    X= coord(:,1)';
                    Y= coord(:,2)';
                    Z= coord(:,3)';

                    plotHexa8(X,Y,Z,'r')
                    set(gca,'Visible','on');
%                     set(gca,'xtick',[],'ytick',[])
                end                               
%%
function plotHexa8(X,Y,Z,c)
plot3([X(1),X(2)],[Y(1),Y(2)],[Z(1),Z(2)],c)
hold on % pour ne pas effacer ce qui a �t� affich� pr�c�demment
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
%
end