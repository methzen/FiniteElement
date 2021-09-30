classdef Plaque
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties % coordonnées nodales et connectivité
        Nodes3D
        Connect3D
        Nodes2D
        Connect2D
        Numerot2D
%         Lcaract
    end
    properties % plaque
        Long
        Larg
        Epaiss
        discretx
        discrety
        discretz
    end
    
    methods
        function obj=Plaque(mLong,mLarg,mEpaiss,mdiscretx,mdiscrety,mdiscretz)
            if nargin==6
                obj.Long=mLong;
                obj.Larg=mLarg;
                obj.Epaiss=mEpaiss;
                obj.discretx=mdiscretx;
                obj.discrety=mdiscrety;
                obj.discretz=mdiscretz;
            [obj.Nodes2D,obj.Connect2D,obj.Numerot2D] = DiscretPlane(obj);
            [obj.Nodes3D,obj.Connect3D] = ExtrudePlane(obj);
            end

        end
         function [N2D,N3D,C3D]=getProperty(obj)              
              N3D=obj.Nodes3D;
              N2D=obj.Numerot2D;
              C3D=obj.Connect3D;
        end
    end
    
end

