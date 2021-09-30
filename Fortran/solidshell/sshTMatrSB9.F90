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
subroutine sshTMatrSB9(nb_node, geom, XI, T)
!
implicit none
!
#include "asterfort/matinv.h"
#include "asterfort/sshJaco.h"
#include "asterfort/assert.h"
!
integer, intent(in) :: nb_node
real(kind=8), intent(in) :: geom(3*nb_node), XI(3)
real(kind=8), intent(out) :: T(6,6)
!
! --------------------------------------------------------------------------------------------------
!
! Solid-shell element - SB9
!
! Compute T matrix (matrix relating the covariant and cartesian frames) at current point
!
! --------------------------------------------------------------------------------------------------
!
! T = [TXI TETA TZETA] = T0 + xi.TXI + eta.TETA + zeta.TZETA
!
! In  nb_node          : number of nodes of element without pinch node(s)
! In  geom             : coordinates of element
! In  XI               : parametric coordinates of current point
! Out T                : matrix relating the covariant and cartesian frames
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: J11, J12, J13, J21, J22, J23, J31, J32, J33, det
    real(kind=8) :: J(3,3), Jm1(3,3)
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(nb_node .eq. 8)
!
    T = 0.d0
!
! - Compute jacobian matrix at current point
!
    call sshJaco(nb_node, geom, XI, J)
!
! - Compute inverse jacobian matrix [J]^-1
!
    call matinv('S', 3, J, Jm1, det)
!
    J11=Jm1(1,1)
    J12=Jm1(1,2)
    J13=Jm1(1,3)
    J21=Jm1(2,1)
    J22=Jm1(2,2)
    J23=Jm1(2,3)
    J31=Jm1(3,1)
    J32=Jm1(3,2)
    J33=Jm1(3,3)
!
    T(1,1)=J11*J11
    T(1,2)=J21*J21
    T(1,3)=J31*J31
    T(1,4)=J11*J21
    T(1,5)=J11*J31
    T(1,6)=J21*J31
!
    T(2,1)=J12*J12
    T(2,2)=J22*J22
    T(2,3)=J32*J32
    T(2,4)=J12*J22
    T(2,5)=J12*J32
    T(2,6)=J22*J32
!
    T(3,1)=J13*J13
    T(3,2)=J23*J23
    T(3,3)=J33*J33
    T(3,4)=J13*J23
    T(3,5)=J13*J33
    T(3,6)=J23*J33
!
    T(4,1)=2.d0*J11*J12
    T(4,2)=2.d0*J21*J22
    T(4,3)=2.d0*J31*J32
    T(4,4)=J11*J22+J21*J12
    T(4,5)=J11*J32+J31*J12
    T(4,6)=J22*J31+J21*J32
!
    T(6,1)=2.d0*J12*J13
    T(6,2)=2.d0*J22*J23
    T(6,3)=2.d0*J32*J33
    T(6,4)=J13*J22+J12*J23
    T(6,5)=J12*J33+J32*J13
    T(6,6)=J22*J33+J32*J23
!
    T(5,1)=2.d0*J11*J13
    T(5,2)=2.d0*J21*J23
    T(5,3)=2.d0*J31*J33
    T(5,4)=J13*J21+J11*J23
    T(5,5)=J13*J31+J11*J33
    T(5,6)=J23*J31+J21*J33
!
end subroutine
