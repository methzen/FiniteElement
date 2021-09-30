%--------------------------------------------------------------------------
% Mise à jour géometrie.
%--------------------------------------------------------------------------
function x = update_geometry(x,eta,displ,freedof)
x(freedof) = x(freedof) + eta*displ;
