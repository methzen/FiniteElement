! --------------------------------------------------------------------
! Copyright (C) 1991 - 2020 - EDF R&D - www.code-aster.org
! This file is part of code_aster.
!
! code_aster is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! code_aster is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with code_aster.  If not, see <http://www.gnu.org/licenses/>.
! --------------------------------------------------------------------
!
subroutine sshGradMatrCartSB9(nb_node  , geom,&
                              nb_dof   , disp     ,&
                              B0       , BZETA    , BZETAZETA,&
                              BXI_     , BETA_    , &
                              BXIZETA_ , BETAZETA_, &
                              det_)
!
implicit none
!
#include "asterf_types.h"
#include "asterfort/sshTMatrDecoSB9.h"
#include "asterfort/sshTMatrSB9.h"
#include "asterfort/sshGradMatrCovaSB9.h"
#include "asterfort/sshJaco.h"
#include "asterfort/matinv.h"
#include "asterfort/assert.h"
!
integer, intent(in) :: nb_node, nb_dof
real(kind=8), intent(in)  :: geom(3*nb_node), disp(nb_dof)
real(kind=8), intent(out) :: B0(6,3*nb_node), BZETA(6,3*nb_node), BZETAZETA(6,3*nb_node)
real(kind=8), optional, intent(out) :: BXI_(6,3*nb_node), BETA_(6,3*nb_node)
real(kind=8), optional, intent(out) :: BXIZETA_(6,3*nb_node), BETAZETA_(6,3*nb_node)
real(kind=8), optional, intent(out) :: det_
!
! --------------------------------------------------------------------------------------------------
!
! Solid-shell element - SB9
!
! Compute the gradient matrix in cartesian frame
!
! --------------------------------------------------------------------------------------------------
!
! In  nb_node          : number of nodes of element without pinch node(s)
! In  geom             : coordinates of element
! In  nb_dof           : number of dof
! In  disp             : displacements of element
! Out B0               : constant part of gradient matrix
! Out BXI              : xi part of gradient matrix
! Out BETA             : eta part of gradient matrix
! Out BZETA            : zeta part of gradient matrix
! Out BXIZETA          : xi x zeta part of gradient matrix
! Out BETAZETA         : eta x zeta part of gradient matrix
! Out BZETAZETA        : zeta x zeta part of gradient matrix
! Out det              : determinant of initial jacobian matrix
!
! NB: from Hughes, the B matrix is already "deviatoric" => B* = B*Dev
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: Bc0(6,24), BcZETA(6,24), BcZETAZETA(6,24), BcXI(6,24), geomCurr(24)
    real(kind=8) :: BcETA(6,24), BcETAZETA(6,24), BcXIZETA(6,24)
    real(kind=8) :: TXI(6,6), TETA(6,6), TZETA(6,6), T(6,6), J0(3,3), Jm1(3,3)
    real(kind=8) :: det
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(nb_node .eq. 8)
    ASSERT(nb_dof .eq. 25)
    B0        = 0.d0
    BZETA     = 0.d0
    BZETAZETA = 0.d0
    det       = 0.d0
    if (present(BXI_)) then
        BXI_      = 0.d0
        BETA_     = 0.d0
        BETAZETA_ = 0.d0
        BXIZETA_  = 0.d0
    endif
!
! - Update configuration
!
    geomCurr(:) = geom(:) + disp(1:3*nb_node)
!
! - Compute jacobian matrix at center of element
!
    call sshJaco(nb_node, geom, (/0.d0, 0.d0, 0.d0/), J0)
!
! - Inverse jacobian matrix
!
    call matinv('S', 3, J0, Jm1, det)
    !WRITE(6,*) 'Jacobien: ',J0,Jm1,det
!
! - Misoriented cell => get absolute !
!
    det = abs(det)
!
! - Compute gradient matrix in covariant basis
!
    call sshGradMatrCovaSB9(nb_node , geomCurr, Bc0       ,&
                            BcXI    , BcETA    , BcZETA    ,&
                            BcXIZETA, BcETAZETA, BcZETAZETA)
!
! - Compute the decomposition of the inverse Jacobian matrix
!
    call sshTMatrDecoSB9(nb_node, geom, TXI, TETA, TZETA)
    !print*,'TXI',TXI
    !print*,'TETA',TETA
    !print*,'TZETA',TZETA
!
! - Compute T matrix (matrix relating the covariant and cartesian frames) at center of element
!
    call sshTMatrSB9(nb_node, geom, (/0.d0,0.d0,0.d0/), T)
    !WRITE(6,*) 'Décomposition jacobien: ',TXI, TETA, TZETA, T
!
! - Decomposition of gradient matrix
!
    B0        = matmul(T,Bc0)
    BZETA     = matmul(T,BcZETA)+matmul(TZETA,Bc0)
    BZETAZETA = matmul(T,BcZETAZETA)+matmul(TZETA,BcZETA)
    if (present(BXI_)) then
        BXI_      = matmul(T,BcXI)+matmul(TXI,Bc0)
        BETA_     = matmul(T,BcETA)+matmul(TETA,Bc0)
        BETAZETA_ = matmul(T,BcETAZETA)+matmul(TETA,BcZETA)+matmul(TZETA,BcETA)
        BXIZETA_  = matmul(T,BcXIZETA)+matmul(TXI,BcZETA)+matmul(TZETA,BcXI)
    endif
!
    if (present(det_)) then
        det_ = det
    endif
!
end subroutine
