%--------------------------------------------------------------------------
% Mise a jour de la force extérieure 
%--------------------------------------------------------------------------
function [Residu,force_ext] = external_force_update(Force_total,...
                                     Residu,force_ext,dlamb)
Residu = Residu - dlamb*Force_total;
force_ext = force_ext + dlamb*Force_total;