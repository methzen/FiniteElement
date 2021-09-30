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
! aslint: disable=W1504
!
subroutine sshGetFaceSB7(Xi, Yi, Zi, Xj, Yj, Zj,&
                         Xm, Ym, Zm, Xl, Yl, Zl,&
                         xS, yS, zS,&
                         xQ, yQ, zQ,&
                         xR, yR, zR,&
                         si, qi, sj, qj,&
                         sl, ql, sm, qm,&
                         faceArea)
!
implicit none
!
real(kind=8), intent(in) :: Xi, Yi, Zi, Xj, Yj, Zj, Xl, Yl, Zl, Xm, Ym, Zm
real(kind=8), intent(out) :: xR, yR, zR
real(kind=8), intent(out) :: xS, yS, zS
real(kind=8), intent(out) :: xQ, yQ, zQ
real(kind=8), intent(out) :: si, qi, sj, qj, sl, ql, sm, qm
real(kind=8), intent(out) :: faceArea
!
! --------------------------------------------------------------------------------------------------
!
! Solid-shell element - SB7
!
! Get face 
!
! --------------------------------------------------------------------------------------------------
!
! Input: quadrangular face from SB7 with nodes i, j, l, m in global coordinates
! Output: quadrangular face with nodes i, j, l, m in (s,r) plane
!
! In  Xi, Yi, Zi       : coordinates of node 1 of quadrangular face
! In  Xj, Yj, Zj       : coordinates of node 2 of quadrangular face
! In  Xl, Yl, Zl       : coordinates of node 3 of quadrangular face
! In  Xm, Ym, Zm       : coordinates of node 4 of quadrangular face
! Out xR, yR, zR       : coordinates of normal to face - R
! Out xS, yS, zS       : coordinates of first tangent of face - S
! Out xQ, yQ, zQ       : coordinates of second tangent of face - Q
! Out si, qi           : local coordinates (in s, q plane) for the node 1 of face
! Out sj, qj           : local coordinates (in s, q plane) for the node 2 of face
! Out sl, ql           : local coordinates (in s, q plane) for the node 3 of face
! Out sm, qm           : local coordinates (in s, q plane) for the node 4 of face
! Out faceArea         : area of face
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: x12, y12, z12, xAB, yAB, zAB
    real(kind=8) :: normR, normQ, normS
!
! --------------------------------------------------------------------------------------------------
!

!
! - Compute horizontal vector in the middle of face
!
    x12 = (Xj+Xm)/2.d0-(Xi+Xl)/2.d0
    y12 = (Yj+Ym)/2.d0-(Yi+Yl)/2.d0
    z12 = (Zj+Zm)/2.d0-(Zi+Zl)/2.d0
!
! - Compute horizontal vector in the middle of face
!
    xAB = (Xl+Xm)/2.d0-(Xi+Xj)/2.d0
    yAB = (Yl+Ym)/2.d0-(Yi+Yj)/2.d0
    zAB = (Zl+Zm)/2.d0-(Zi+Zj)/2.d0
!
! - Normal to face r = 12 ^ AB (normalized)
!
    xR = y12*zAB-z12*yAB
    yR = z12*xAB-x12*zAB
    zR = x12*yAB-y12*xAB
    normR = sqrt(xR*xR+yR*yR+zR*zR)
    xR = xR/normR
    yR = yR/normR
    zR = zR/normR
!
! - First tangent to face s = 12 (normalized)
!
    normS = sqrt(x12*x12+y12*y12+z12*z12)
    xS = x12/normS
    yS = y12/normS
    zS = z12/normS
!
! - Second tangent to face q = r ^ s (normalized)
!
    xQ = yR*zS-zR*yS;
    yQ = zR*xS-xR*zS;
    zQ = xR*yS-yR*xS;
    normQ  = sqrt(xQ*xQ+yQ*yQ+zQ*zQ);
    xQ = xQ/normQ
    yQ = yQ/normQ
    zQ = zQ/normQ
!
! - Local coordinates of face (in s,q plane)
!
    si = xS*Xi+yS*Yi+zS*Zi
    qi = xQ*Xi+yQ*Yi+zQ*Zi
    sj = xS*Xj+yS*Yj+zS*Zj
    qj = xQ*Xj+yQ*Yj+zQ*Zj
    sm = xS*Xm+yS*Ym+zS*Zm
    qm = xQ*Xm+yQ*Ym+zQ*Zm
    sl = xS*Xl+yS*Yl+zS*Zl
    ql = xQ*Xl+yQ*Yl+zQ*Zl
!
! - Area of face
!
    faceArea = abs(((sm-si)*(ql-qj)+(qm-qi)*(sj-sl)))/2.d0
!
end subroutine
