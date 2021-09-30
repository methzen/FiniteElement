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
! aslint: disable=W1306
!
subroutine sshNLMatr(nbsig   , nb_node, nb_dof,&
                     jacobian, B      , dsidep,&
                     matuu   , kGeom_)
!
implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/sshBTSB.h"
!
integer, intent(in) :: nb_node, nb_dof, nbsig
real(kind=8), intent(in) :: jacobian
real(kind=8), intent(in) :: dsidep(nbsig, nbsig), B(nbsig, nb_dof)
real(kind=8), intent(inout) :: matuu(nb_dof*nb_dof)
real(kind=8), optional, intent(in) :: kGeom_(nb_dof, nb_dof)
!
! --------------------------------------------------------------------------------------------------
!
! Solid-shell element
!
! Compute matrix
!
! --------------------------------------------------------------------------------------------------
!
! In  nbsig            : number of components o the stress tensor
! In  nb_node          : number of nodes of element without pinch node(s)
! In  nb_dof           : number of dof
! In  jacobian         : product of determinant and gauss weight
! In  B                : gradient B matrix (reduced integration)
! In  dsidep           : jacobian matrix of behaviour (dSigm/dEpsi)
! In  kGeom            : geometric part of matrix
! IO  matuu            : elastic matrix
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8)  :: btdb(nb_dof, nb_dof)
    integer :: i, j, k
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(nbsig .eq. 6)
    ASSERT(nb_node .eq. 8 .or. nb_node .eq. 6)
    ASSERT(nb_dof .eq. 25 .or. nb_dof .eq. 19)
!
! - Compute material part
!
    call sshBTSB(dsidep, 6, nb_dof, B, btdb)
    !WRITE(6,*) 'Matrice matérielle: ',sum(btdb)
!
! - Add geometric part
!
    if (present(kGeom_)) then
        do i = 1, 3*nb_node
            do j = 1, 3*nb_node
                btdb(i,j) = btdb(i,j) + kGeom_(i,j)
            end do
        end do
        !WRITE(6,*) 'Matrice géométrique: ',sum(kGeom_)
    endif
!
! - Compute tangent matrix
!
    k = 0
    do i = 1, nb_dof
        do j = 1, i
            k = k + 1
            matuu(k) = matuu(k) + jacobian*btdb(i,j)
        enddo
    enddo
    !WRITE(6,*) 'Jacobien: ',jacobian
    !WRITE(6,*) 'Matrice B: ',sum(B)
    !WRITE(6,*) 'Matrice avant stabilisation: ',sum(matuu(1:nb_dof*nb_dof))
!
end subroutine
