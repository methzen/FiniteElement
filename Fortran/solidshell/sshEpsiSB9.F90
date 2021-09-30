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
subroutine sshEpsiSB9(nb_node , nb_dof,&
                      zeta    , disp  ,&
                      B0      , BZETA , BZETAZETA,&
                      epsi)
!
use Behaviour_type
!
implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
!
integer, intent(in) :: nb_node, nb_dof
real(kind=8), intent(in) :: disp(nb_dof), zeta
real(kind=8), intent(in) :: B0(6, 3*nb_node), BZETA(6, 3*nb_node), BZETAZETA(6, 3*nb_node)
real(kind=8), intent(out):: epsi(6)
!
! --------------------------------------------------------------------------------------------------
!
! Solid-shell element - SB9
!
! Compute small strains in cartesian frame at current point
!
! --------------------------------------------------------------------------------------------------
!
!
! In  l_large          : flag for large strains
! In  nb_node          : number of nodes of element without pinch node(s)
! In  nb_dof           : number of dof
! In  geomInit         : initial coordinates of element
! In  disp             : displacements of element
! In  zeta             : out-of-plane parametric component
! In  B0               : constant part of gradient matrix
! In  BZETA            : zeta part of gradient matrix
! In  BZETAZETA        : zeta x zeta part of gradient matrix
! Out epsi             : small strains
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(nb_node .eq. 8)
    ASSERT(nb_dof .eq. 25)
    epsi = 0.d0
    epsi = matmul((B0+BZETA*zeta+BZETAZETA*zeta*zeta), disp(1:3*nb_node))
!
end subroutine
