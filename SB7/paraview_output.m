%-----------------------------------------------------------------------
function NOPRVW=paraview_output(E,Xcoord,DEP,NOPRVW)
%-----------------------------------------------------------------------
%
%     name   --> main part of the output file name, used for paraview files!
%
%-----------------------------------------------------------------------

%
%% -----Based upon the "incout" write out the mesh and results for paraview plotting
% 
ndime=3;
%Nd=ChDeplacement(obj);
        %N=full(obj.Nodes3D);
        %E=obj.Connect3D;
%
E=E-1;  % on fait commencer les valeurs � 0 plut�t que 0
Nd=DEP;

%Xcoord=N(:,2:4);
    if isempty(dir('OutputParaviewFolder'))
        mkdir('./','OutputParaviewFolder')
    end
    npoin=length(Xcoord);   %  number of nodes
     [nelem,NbrElemNode]=size(E);    %  number of elements
    if isnumeric(NOPRVW)
        exten=sprintf('%3.3i',NOPRVW);  
    else
        exten=NOPRVW;
    end
    PARname=['./','OutputParaviewFolder/','paraMatlab',exten,'.vtk'];
    PARfile=fopen(PARname,'w') ;
    
    
    
    fprintf(PARfile, '# vtk DataFile Version 3.0 \r\n');
    fprintf(PARfile,'vtk output \r\n');
    fprintf(PARfile, 'ASCII \r\n');
    fprintf(PARfile,'DATASET UNSTRUCTURED_GRID \r\n');
%%	ecriture des coorfornn�es generalis�	
    fprintf(PARfile,strcat(sprintf('POINTS    %10i     ',npoin),' float\r\n'));
	%fprintf(PARfile,'0 0 0 \r\n');

    for J=1:npoin
        if (ndime==3)
            fprintf(PARfile,strcat(sprintf('%12.4E  ',Xcoord(J,1:3)),'\r\n'));
        else
            fprintf(PARfile,strcat(sprintf('%12.4E  ',Xcoord(J,1:2)),'       0 \r\n'));
        end
    end
    
    
%%	ecriture de la tables des connectivit�s
% calcul du nombre de valeur qui vont �tre introduite 

NbrOfData=nelem*(NbrElemNode+1);

fprintf(PARfile,'CELLS  %10i %10i \r\n',nelem,NbrOfData );

for J=1:nelem  
fprintf(PARfile,' %10i  %10i  %10i  %10i  %10i  %10i  %10i  %10i  %10i \r\n', [8 E(J,:)]);
end

%%	ecriture du type d'�l�ments
    
fprintf(PARfile,strcat(sprintf('CELL_TYPES     %s',int2str(nelem)),'\r\n'));
for J=1:nelem

         fprintf(PARfile,' %i \r\n',12);
       
      
end

Dep=Nd;
%%	le champ de d�placement	
    fprintf(PARfile,strcat(sprintf('POINT_DATA    %10i     ',npoin),'  \r\n'));
    fprintf(PARfile,'VECTORS deplacement   float64  \r\n');
	%fprintf(PARfile,'0 0 0 \r\n');
    for J=1:npoin
        k=(J-1)*3+1;
        if (ndime==3)
            UU=[Dep(k),Dep(k+1),Dep(k+2)];
            fprintf(PARfile,strcat(sprintf('%12.4E  ',UU),'\r\n'));
        else
            fprintf(PARfile,strcat(sprintf('%12.4E  ',Dep(J,1:2)),'       0 \r\n'));
        end
    end
%%	le champ des d�formations	  

%%
fclose(PARfile);
if isnumeric(NOPRVW)
    NOPRVW=NOPRVW+1;
elseif strcmpi(NOPRVW,'Initial')
    NOPRVW=0;
end

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end 
