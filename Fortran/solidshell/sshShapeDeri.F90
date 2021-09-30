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
subroutine sshShapeDeri(nb_node, XI, dN_dXsi, dN_dEta, dN_dZeta)
!
implicit none
!
#include "asterfort/assert.h"
!
integer, intent(in) :: nb_node
real(kind=8), intent(in)  :: XI(3)
real(kind=8), intent(out) :: dN_dXsi(nb_node), dN_dEta(nb_node), dN_dZeta(nb_node)
!
! --------------------------------------------------------------------------------------------------
!
! Solid-shell element
!
! Compute the decomposed derivatives of shape functions at current point
!
! --------------------------------------------------------------------------------------------------
!
! In  nb_node          : number of nodes of element without pinch node(s)
! In  XI               : parametric coordinates of current point
! In  dN_dXsi          : derivatives of shapes function by XI
! In  dN_dEta          : derivatives of shapes function by ETA
! In  dN_dZeta         : derivatives of shapes function by DZETA
!
! --------------------------------------------------------------------------------------------------
!
    integer :: i
    real(kind=8) :: g1(8), g2(8), g3(8)
    real(kind=8) :: h1(8), h2(8), h3(8), h4(8)
    real(kind=8), parameter :: un = 1.d0, demi = 0.5d0, zero = 0.d0
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(nb_node .eq. 8 .or. nb_node .eq. 6)
!
    dN_dXsi  = 0.d0
    dN_dEta  = 0.d0
    dN_dZeta = 0.d0

    if (nb_node .eq. 8) then 
!
! - Vectors for decomposition
!
        g1=(/-1.d0,  1.d0,  1.d0, -1.d0, -1.d0,  1.d0, 1.d0, -1.d0/)/8.d0
        g2=(/-1.d0, -1.d0,  1.d0,  1.d0, -1.d0, -1.d0, 1.d0,  1.d0/)/8.d0
        g3=(/-1.d0, -1.d0, -1.d0, -1.d0,  1.d0,  1.d0, 1.d0,  1.d0/)/8.d0
        h1=(/ 1.d0, -1.d0,  1.d0, -1.d0,  1.d0, -1.d0, 1.d0, -1.d0/)/8.d0
        h2=(/ 1.d0,  1.d0, -1.d0, -1.d0, -1.d0, -1.d0, 1.d0,  1.d0/)/8.d0
        h3=(/ 1.d0, -1.d0, -1.d0,  1.d0, -1.d0,  1.d0, 1.d0, -1.d0/)/8.d0
        h4=(/-1.d0,  1.d0, -1.d0,  1.d0,  1.d0, -1.d0, 1.d0, -1.d0/)/8.d0
!
        do i = 1, nb_node
            dN_dXsi(i)  = g1(i)+XI(2)*h1(i)+XI(3)*h3(i)+XI(2)*XI(3)*h4(i)
            dN_dEta(i)  = g2(i)+XI(1)*h1(i)+XI(3)*h2(i)+XI(1)*XI(3)*h4(i)
            dN_dZeta(i) = g3(i)+XI(2)*h2(i)+XI(1)*h3(i)+XI(2)*XI(1)*h4(i)
        enddo

    else
        dN_dXsi(1)=demi*(-1.d0)*(un-XI(3))
        dN_dXsi(2)=demi*(+1.d0)*(un-XI(3))
        dN_dXsi(3)=zero
        dN_dXsi(4)=demi*(-1.d0)*(un+XI(3))
        dN_dXsi(5)=demi*(+1.d0)*(un+XI(3))
        dN_dXsi(6)=zero
!
        dN_dEta(1)=demi*(-1.d0)*(un-XI(3))
        dN_dEta(2)=zero
        dN_dEta(3)=demi*(+1.d0)*(un-XI(3))
        dN_dEta(4)=demi*(-1.d0)*(un+XI(3))
        dN_dEta(5)=zero
        dN_dEta(6)=demi*(+1.d0)*(un+XI(3))
!
        dN_dZeta(1)=demi*(un-XI(1)-XI(2))*(-1.d0)
        dN_dZeta(2)=demi*XI(1)*(-1.d0)
        dN_dZeta(3)=demi*XI(2)*(-1.d0)
        dN_dZeta(4)=demi*(un-XI(1)-XI(2))*(+1.d0)
        dN_dZeta(5)=demi*XI(1)*(+1.d0)
        dN_dZeta(6)=demi*XI(2)*(+1.d0)
    endif
!
end subroutine
