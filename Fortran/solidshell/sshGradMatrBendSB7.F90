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
subroutine sshGradMatrBendSB7(nb_node , geom, R,&
                              triaN1x , triaN1y,&
                              triaN2x , triaN2y,&
                              triaN3x , triaN3y,&
                              triaArea, Bb)
!
implicit none
!
#include "asterfort/assert.h"
#include "asterfort/dkt6GradMatrBend.h"
#include "asterfort/sshRotaDispSB7.h"
!
integer, intent(in) :: nb_node
real(kind=8), intent(in) :: geom(3*nb_node), R(3,3)
real(kind=8), intent(in) :: triaN1x, triaN1y, triaN2x, triaN2y, triaN3x, triaN3y
real(kind=8), intent(in) :: triaArea
real(kind=8), intent(out) :: Bb(3,18)
!
! --------------------------------------------------------------------------------------------------
!
! Solid-shell element - SB7
!
! Compute bending terms of the gradient B matrix
!
! --------------------------------------------------------------------------------------------------
!
! In  nb_node          : number of nodes of element without pinch node(s)
! In  geom             : coordinates of element
! In  R                : rotation matrix for frame of mid_edge triangle
! In  triaN1x, triaN1y : local coordinates of mid_edge triangle, first vertex
! In  triaN2x, triaN2y : local coordinates of mid_edge triangle ,second vertex
! In  triaN3x, triaN3y : local coordinates of mid_edge triangle, third vertex
! In  triaArea           : area of mid_edge triangle
! Out Bb               : bending part of the gradient B matrix
!
! --------------------------------------------------------------------------------------------------
!
    integer :: p
    real(kind=8) :: Bbw(3,18),Bt(3,3),Bw(3,3)
    real(kind=8) :: Bbt(3,18), Tt(3, 18)
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(nb_node .eq. 6)
    Bb = 0.d0
!
! - Compute THETA and W terms (rotation/displacement) from DKT6 element theory
!
    !WRITE(6,*) 'DKT6-In',triaN1x , triaN1y, triaN2x,&
    !              triaN2y , triaN3x, triaN3y
    call dkt6GradMatrBend(triaN1x , triaN1y, triaN2x,&
                          triaN2y , triaN3x, triaN3y,&
                          triaArea, Bt  , Bw)
    !WRITE(6,*) 'DKT6',triaArea,sum(Bt),sum(Bw)
!
! - Compute bending terms from rotations
!
    call sshRotaDispSB7(geom, Tt)
    !WRITE(6,*) 'Bending Tt',sum(Tt)
    Bbt = matmul(Bt, Tt)
    !WRITE(6,*) 'Bending Bb',sum(Bbt)
!
! - Compute bending terms from displacements
!
    Bbw = 0.d0
    do p = 1, 3
        Bbw(p,1)  = 0.5d0*Bw(p,1)*R(3,1)
        Bbw(p,2)  = 0.5d0*Bw(p,1)*R(3,2)
        Bbw(p,3)  = 0.5d0*Bw(p,1)*R(3,3)
        Bbw(p,10) = Bbw(p,1)
        Bbw(p,11) = Bbw(p,2)
        Bbw(p,12) = Bbw(p,3)
        Bbw(p,4)  = 0.5d0*Bw(p,2)*R(3,1)
        Bbw(p,5)  = 0.5d0*Bw(p,2)*R(3,2)
        Bbw(p,6)  = 0.5d0*Bw(p,2)*R(3,3)
        Bbw(p,13) = Bbw(p,4)
        Bbw(p,14) = Bbw(p,5)
        Bbw(p,15) = Bbw(p,6)
        Bbw(p,7)  = 0.5d0*Bw(p,3)*R(3,1)
        Bbw(p,8)  = 0.5d0*Bw(p,3)*R(3,2)
        Bbw(p,9)  = 0.5d0*Bw(p,3)*R(3,3)
        Bbw(p,16) = Bbw(p,7)
        Bbw(p,17) = Bbw(p,8)
        Bbw(p,18) = Bbw(p,9)
    enddo
    !WRITE(6,*) 'Bending Bbw',sum(Bbw)
!
! - Compute bending
!
    Bb = Bbt + Bbw
!
end subroutine
