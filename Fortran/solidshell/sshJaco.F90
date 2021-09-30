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
! aslint: disable=W1306
!
subroutine sshJaco(nb_node, geom, XI, Jac)
!
implicit none
!
#include "asterfort/assert.h"
#include "asterfort/sshShapeDeri.h"
!
integer, intent(in) :: nb_node
real(kind=8), intent(in)  :: geom(3*nb_node), XI(3)
real(kind=8), intent(out) :: Jac(3,3)
!
! --------------------------------------------------------------------------------------------------
!
! Solid-shell element
!
! Compute jacobian matrix at current point
!
! --------------------------------------------------------------------------------------------------
!
! In  nb_node          : number of nodes of element without pinch node(s)
! In  geom             : coordinates of element
! In  XI               : parametric coordinates of current point
! Out Jac              : jacobian matrix
!
! --------------------------------------------------------------------------------------------------
!
    integer :: i
    real(kind=8) :: dN_dXsi(nb_node), dN_dEta(nb_node), dN_dZeta(nb_node)
    real(kind=8) :: dX_dXsi, dX_dEta, dX_dZeta
    real(kind=8) :: dZ_dXsi, dZ_dEta, dZ_dZeta
    real(kind=8) :: dY_dXsi, dY_dEta, dY_dZeta
    real(kind=8) :: X(nb_node), Y(nb_node), Z(nb_node)
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(nb_node .eq. 8 .or. nb_node .eq. 6)
!
    do i = 1, nb_node
        X(i) = geom(3*(i-1)+1)
        Y(i) = geom(3*(i-1)+2)
        Z(i) = geom(3*(i-1)+3)
    enddo
!
! - Derivatives of shape functions
!
    call sshShapeDeri(nb_node, XI, dN_dXsi, dN_dEta, dN_dZeta)
!
! - Compute jacobian
!
    dX_dXsi  = sum(X*dN_dXsi)
    dX_dEta  = sum(X*dN_dEta)
    dX_dZeta = sum(X*dN_dZeta)
!
    dY_dXsi  = sum(Y*dN_dXsi)
    dY_dEta  = sum(Y*dN_dEta)
    dY_dZeta = sum(Y*dN_dZeta)
!
    dZ_dXsi  = sum(Z*dN_dXsi)
    dZ_dEta  = sum(Z*dN_dEta)
    dZ_dZeta = sum(Z*dN_dZeta)
!
    Jac(1,1) = dX_dXsi
    Jac(1,2) = dY_dXsi
    Jac(1,3) = dZ_dXsi
    Jac(2,1) = dX_dEta
    Jac(2,2) = dY_dEta
    Jac(2,3) = dZ_dEta
    Jac(3,1) = dX_dZeta
    Jac(3,2) = dY_dZeta
    Jac(3,3) = dZ_dZeta
    Jac = transpose(Jac)
!
end
