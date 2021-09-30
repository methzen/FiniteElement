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
subroutine sshGradMatrCovaSB7(nb_node, geomLocal, XI, bx, by, bz)
!
implicit none
!
#include "asterfort/assert.h"
#include "asterfort/sshJaco.h"
#include "asterfort/matinv.h"
!
integer, intent(in) :: nb_node
real(kind=8), intent(in) :: geomLocal(3*nb_node), XI(3)
real(kind=8), intent(out) :: bx(6), by(6), bz(6)
!
! --------------------------------------------------------------------------------------------------
!
! Solid-shell element - SB7
!
! Compute gradient matrix in covariant basis
!
! --------------------------------------------------------------------------------------------------
!
! In  nb_node          : number of nodes of element without pinch node(s)
! In  geomLocal        : coordinates of element in the frame of mid_edge triangle
! In  XI               : parametric coordinates where evaluate gradient matrix
! Out bx               : gradient matrix in local covariant basis - X coordinate
! Out by               : gradient matrix in local covariant basis - Y coordinate
! Out bz               : gradient matrix in local covariant basis - Z coordinate
!
! --------------------------------------------------------------------------------------------------
!
    integer :: j
    real(kind=8) :: Bksi(3,6)
    real(kind=8), parameter :: uns2 = 1.d0/2.d0, uns6 = 1.d0/6.d0
    real(kind=8) :: J1s3(3,3), invJ1s3(3,3), detJ
!
! --------------------------------------------------------------------------------------------------
!
    bx = 0.d0
    by = 0.d0
    bz = 0.d0
!
! - Local jacobian 
!
    call sshJaco(nb_node, geomLocal, XI, J1s3)
    call matinv('S', 3, transpose(J1s3), invJ1s3, detJ)
!
! - <N,xi> <N,eta> <N,zeta>--|(1/3,1/3,zeta)
!
    Bksi = 0.d0
    Bksi(1,1) = -uns2+uns2*XI(3)
    Bksi(1,2) =  uns2-uns2*XI(3)
    Bksi(1,3) = 0.d0
    Bksi(1,4) = -uns2-uns2*XI(3)
    Bksi(1,5) = uns2+uns2*XI(3)
    Bksi(1,6) = 0.d0
!
    Bksi(2,1) = -uns2+uns2*XI(3)
    Bksi(2,2) = 0.d0
    Bksi(2,3) = uns2-uns2*XI(3)
    Bksi(2,4) = -uns2-uns2*XI(3)
    Bksi(2,5) = 0.d0
    Bksi(2,6) = uns2+uns2*XI(3)
!
    Bksi(3,1) = -uns6
    Bksi(3,2) = -uns6
    Bksi(3,3) = -uns6
    Bksi(3,4) = uns6
    Bksi(3,5) = uns6
    Bksi(3,6) = uns6
!
! - line 1:<N,x>  line 2:<N,y>  line 3:<N,z> - en xi=1/3, eta=1/3, zeta=0
!
    do j = 1, 6
        bx(j) = invJ1s3(1,1)*Bksi(1,j)+invJ1s3(1,2)*Bksi(2,j)+invJ1s3(1,3)*Bksi(3,j)
        by(j) = invJ1s3(2,1)*Bksi(1,j)+invJ1s3(2,2)*Bksi(2,j)+invJ1s3(2,3)*Bksi(3,j)
        bz(j) = invJ1s3(3,1)*Bksi(1,j)+invJ1s3(3,2)*Bksi(2,j)+invJ1s3(3,3)*Bksi(3,j)
    end do
!
end subroutine
