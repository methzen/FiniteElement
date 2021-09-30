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
subroutine sshGradMatrEAS(nb_node, geom, zeta, B9)
!
implicit none
!
#include "asterfort/sshTMatrSB9.h"
!
integer, intent(in) :: nb_node
real(kind=8), intent(in) :: geom(3*nb_node), zeta
real(kind=8), intent(out) :: B9(6)
!
! --------------------------------------------------------------------------------------------------
!
! Solid-shell element - SB9
!
! Compute component gradient for EAS effect
!
! --------------------------------------------------------------------------------------------------
!
! In  nb_node          : number of nodes of element without pinch node(s)
! In  geom             : coordinates of element
! In  zeta             : out-of-plane parametric component
! Out B9               : gradient matrix for EAS effect
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: Bc9(6), T(6,6)
!
! --------------------------------------------------------------------------------------------------
!
    B9  = 0.d0
    call sshTMatrSB9(nb_node, geom, (/0.d0,0.d0,0.d0/), T)
    Bc9 = (/0.d0, 0.d0,-2.d0*zeta, 0.d0, 0.d0, 0.d0/)
    B9  = matmul(T, Bc9)
!
end subroutine
