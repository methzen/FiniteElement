%--------------------------------------------------------------------------
% Check for equilibrium convergence.
%--------------------------------------------------------------------------
function rnorm = check_residual_norm(RESIDU,lamb,FORCEXT,freedof,FIXEDDOF)
%--------------------------------------------------------------------------
% Obtain the reactions.
%--------------------------------------------------------------------------
REACT=RESIDU(FIXEDDOF)+FORCEXT(FIXEDDOF);
%--------------------------------------------------------------------------
% Evaluates the residual norm.
%--------------------------------------------------------------------------
rnorm = dot(RESIDU(freedof),RESIDU(freedof));
fnorm = dot(FORCEXT(freedof),FORCEXT(freedof));
fnorm = fnorm*lamb^2;
enorm =  dot(REACT,REACT);
rnorm =  sqrt(rnorm/(fnorm + enorm));
%--------------------------------------------------------------------------      

