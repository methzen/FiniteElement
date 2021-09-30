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
subroutine sshNLVect(nbsig   , nb_node, nb_dof,&
                     jacobian, B      , sigm  ,&
                     vectu)
!
implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
!
integer, intent(in) :: nb_node, nbsig, nb_dof
real(kind=8), intent(in) :: jacobian
real(kind=8), intent(in) :: B(nbsig, nb_dof)
real(kind=8), intent(in) :: sigm(nbsig)
real(kind=8), intent(inout) :: vectu(nb_dof)
!
! --------------------------------------------------------------------------------------------------
!
! Solid-shell element
!
! Compute internal force at current Gauss point
!
! --------------------------------------------------------------------------------------------------
!
! In  nbsig            : number of components of the stress tensor
! In  nb_node          : number of nodes of element without pinch node(s)
! In  nb_dof           : number of dof
! In  jacobian         : product of determinant and gauss weight
! In  B                : gradient B matrix (reduced integration)
! In  sigm             : Cauchy stresses
! IO  matuu            : elastic matrix
!
! --------------------------------------------------------------------------------------------------
!
    integer :: i, isig
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(nbsig .eq. 6)
    ASSERT(nb_node .eq. 8 .or. nb_node .eq. 6)
    ASSERT(nb_dof .eq. 25 .or. nb_dof .eq. 19)
!
! - Compute internal force
!
    do i = 1, nb_dof
        do isig = 1, 6
            vectu(i) = vectu(i) + B(isig,i)*sigm(isig)*jacobian
        enddo
    enddo
    !WRITE(6,*) 'Vecteur avant stabilisation: ',sum(vectu(1:nb_dof))
!
end subroutine
