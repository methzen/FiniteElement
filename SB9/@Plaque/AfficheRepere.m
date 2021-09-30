function AfficheRepere(obj,x,y,z,scale)
%UNTITLED2 Summary of this function goes here
%   Affiche le repère dans le maillage
hold on
% axe X

    X=quiver3(x,y,z,scale,0,0);
    X.Color='g';
    X.LineWidth = 3;
    X.AutoScale = 'on';
    text(x+scale+0.001,y+0.001,z+0.001,{'X'})
% axe Y
    Y=quiver3(x,y,z,0,scale,0);
    Y.Color='y';
    Y.LineWidth = 3;
    Y.AutoScale = 'on';
    text(x+0.001,y+scale+0.001,z+0.001,{'Y'})
% axe Z
    Z=quiver3(x,y,z,0,0,scale);
    Z.Color='r';
    Z.LineWidth = 3;
    Z.AutoScale = 'on';
    text(x+0.001,y+0.001,z+scale+0.001,{'Z'})
end

