function PrintMesh(Nu,XYZ, Ec,scale)
                [m,n]=size(Ec);
                coord=zeros(6,3);
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
                    XYZi=[X' Y' Z'];
%                     if n==6 
                    plotHexa6(XYZi,'b')
%                     else
%                     plotHexa8(XYZi,'b')
%                     end
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
                    XYZf=[X' Y' Z'];
                    if n==6 
                    plotHexa6(XYZf,'r')
                    else
                    plotHexa8(XYZf,'r')
                    end
                    set(gca,'Visible','on');
%                     set(gca,'xtick',[],'ytick',[])
                end                               
%%
end