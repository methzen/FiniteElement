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
subroutine sshNLStabSB7(lMatr   , lVect ,&
                        nb_node , nb_dof,&
                        jacobian, Ueff  ,&
                        dispIncr, BStab ,&
                        matuu   , vectu)
!
implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/sshBTSB.h"
!
aster_logical, intent(in) :: lMatr, lVect
integer, intent(in) :: nb_node, nb_dof
real(kind=8), intent(in) :: dispIncr(nb_dof)
real(kind=8), intent(in)  :: Ueff, jacobian
real(kind=8), intent(in) :: BStab(2, 18)
real(kind=8), intent(inout) :: vectu(*), matuu(*)
!
! --------------------------------------------------------------------------------------------------
!
! Solid-shell element - SB7
!
! Compute stabilization (twist mode) for matrix and internal force
!
! --------------------------------------------------------------------------------------------------
!
! In  lMatr            : flag to compute matrix
! In  lVect            : flag to compute vector
! In  BStab            : pinch part stabilization of the gradient B matrix
! In  dispIncr         : increment of displacement
! In  jacobian         : product of determinant and gauss weight
! In  Ueff             : effective shear modulus
! IO  vectu            : linear internal force vector
! IO  matuu            : elastic matrix
!
! --------------------------------------------------------------------------------------------------
!
    integer :: i, j, k
    real(kind=8) :: tBsDBs(18, 18)
    real(kind=8), parameter :: eye(2,2) = reshape((/1.d0, 0.d0,&
                                                    0.d0, 1.d0/),(/2,2/))
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(nb_node .eq. 6)
    ASSERT(nb_dof .eq. 19)
!
    call sshBTSB(eye*Ueff, 2, 3*nb_node, Bstab, tBsDBs)
    !WRITE(6,*) 'Matrice de stabilité: ',jacobian,sum(tBsDBs)
    if (lMatr) then
        k = 0
        do i = 1, 3*nb_node
            do j = 1, i
                k = k + 1
                matuu(k) = matuu(k) + jacobian*tBsDBs(i,j)
            enddo
        enddo
    endif
    if (lVect) then
        do i = 1, 3*nb_node
            do j = 1, 3*nb_node
                vectu(i) = vectu(i) + tBsDBs(i,j)*dispIncr(j)*jacobian
            enddo
        enddo
    endif
!
end subroutine
