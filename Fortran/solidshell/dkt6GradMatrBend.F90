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
subroutine dkt6GradMatrBend(x1  , y1, x2,&
                            y2  , x3, y3,&
                            area, Bt, Bw)
!
implicit none
!
real(kind=8), intent(in) :: x1, y1, x2, y2, x3, y3, area
real(kind=8), intent(out) :: Bt(3,3), Bw(3,3)
!
! --------------------------------------------------------------------------------------------------
!
! Solid-shell element - SB7 (but from DKT6 theory)
!
! Compute bending terms of the gradient B matrix - Rotation and displacement terms (from DKT6)
!
! --------------------------------------------------------------------------------------------------
!
! In  x1, y1           : coordinates of node 1 for DKT6
! In  x2, y2           : coordinates of node 2 for DKT6
! In  x3, y3           : coordinates of node 3 for DKT6
! In  area             : area of DKT6
! Out Bt               : B matrix for bending - Rotation terms
! Out Bw               : B matrix for bending - Displacements terms
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: x21, y21, L21
    real(kind=8) :: x32, y32, L32
    real(kind=8) :: x13, y13, L13
    real(kind=8) :: c4, s4, c5, s5, c6, s6
!
! --------------------------------------------------------------------------------------------------
!
    Bt = 0.d0
    Bw = 0.d0
!
! - Geometric quantities
!
    x21=x2-x1
    y21=y2-y1
    L21=sqrt(x21*x21+y21*y21)
    c4=x21/L21
    s4=y21/L21
    x32=x3-x2
    y32=y3-y2
    L32=sqrt(x32*x32+y32*y32)
    c5=x32/L32
    s5=y32/L32
    x13=x1-x3
    y13=y1-y3
    L13=sqrt(x13*x13+y13*y13)
    c6=x13/L13
    s6=y13/L13
!
! - B matrix - Displacement terms
!
    Bw(1,1)=(s4*c4-s6*c6)/area
    Bw(1,2)=(s5*c5-s4*c4)/area
    Bw(1,3)=(s6*c6-s5*c5)/area
    Bw(2,1)=-Bw(1,1)
    Bw(2,2)=-Bw(1,2)
    Bw(2,3)=-Bw(1,3)
    Bw(3,1)=(s4*s4-c4*c4+c6*c6-s6*s6)/area
    Bw(3,2)=(s5*s5-c5*c5+c4*c4-s4*s4)/area
    Bw(3,3)=(s6*s6-c6*c6+c5*c5-s5*s5)/area
!
! - B matrix - Rotation terms
!
    Bt(1,1)=(y21)*s4/area
    Bt(1,2)=(y32)*s5/area
    Bt(1,3)=(y13)*s6/area
    Bt(2,1)=(x21)*c4/area
    Bt(2,2)=(x32)*c5/area
    Bt(2,3)=(x13)*c6/area
    Bt(3,1)=((-y21)*c4+(-x21)*s4)/area
    Bt(3,2)=((-y32)*c5+(-x32)*s5)/area
    Bt(3,3)=((-y13)*c6+(-x13)*s6)/area
!
end subroutine
