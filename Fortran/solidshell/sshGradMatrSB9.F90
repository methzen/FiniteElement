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
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with code_aster.  If not, see <http://www.gnu.org/licenses/>.
! --------------------------------------------------------------------
!
subroutine sshGradMatrSB9(nb_node, nb_dof,&
                          zeta   , geom  ,&
                          B0     , BZETA , BZETAZETA,&
                          B      , B9_)
!
implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/sshGradMatrEAS.h"
!
integer, intent(in) :: nb_node, nb_dof
real(kind=8), intent(in) :: geom(3*nb_node), zeta
real(kind=8), intent(in) :: B0(6, 3*nb_node), BZETA(6, 3*nb_node), BZETAZETA(6, 3*nb_node)
real(kind=8), intent(out) :: B(6, nb_dof)
real(kind=8), optional, intent(out) :: B9_(6)
!
! --------------------------------------------------------------------------------------------------
!
! Solid-shell element - SB9
!
! Compute B matrix at current Gauss point
!
! --------------------------------------------------------------------------------------------------
!
! In  nb_node          : number of nodes of element without pinch node(s)
! In  nb_dof           : number of dof
! In  zeta             : out-of-plane parametric component
! In  geom             : coordinates of element
! In  B0               : constant part of gradient matrix
! In  BZETA            : zeta part of gradient matrix
! In  BZETAZETA        : zeta x zeta part of gradient matrix
! Out B                : gradient B matrix
! Out B9               : gradient matrix for EAS effect
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: B9(6)
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(nb_node .eq. 8)
    ASSERT(nb_dof .eq. 25)
    B  = 0.d0
    B9 = 0.d0
!
! - Compute component gradient for EAS effect
!
    call sshGradMatrEAS(nb_node, geom, zeta, B9)
!
! - Construct B matrix at current Gauss point
!
    B(1:6, 1:24) = B0+zeta*BZETA+zeta*zeta*BZETAZETA
    B(1:6, 25)   = B9
!
    if (present(B9_)) then
        B9_ = B9
    endif
!
end subroutine
