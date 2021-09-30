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
subroutine sshRotaDispSB7(geom, T)
!
implicit none
!
#include "asterfort/sshGetFaceSB7.h"
!
real(kind=8), intent(in)  :: geom(18)
real(kind=8), intent(out) :: T(3, 18)
!
! --------------------------------------------------------------------------------------------------
!
! Solid-shell element - SB7 (but from DKT6 theory)
!
! Get rotation from displacements
!
! --------------------------------------------------------------------------------------------------
!
! In  geom             : coordinates of SB7 element
! Out T                : displacements of nodes from rotations
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: Xi, Yi, Zi, Xj, Yj, Zj, Xk, Yk, Zk, Xl, Yl, Zl, Xm, Ym, Zm, Xn, Yn, Zn
    real(kind=8) :: xR, yR, zR
    real(kind=8) :: xS, yS, zS
    real(kind=8) :: xQ, yQ, zQ
    real(kind=8) :: si, qi, sj, qj, sk, qk, sl, ql, sm, qm, sn, qn
    real(kind=8) :: faceArea
    real(kind=8) :: Pf, Qf
!
! --------------------------------------------------------------------------------------------------
!
    T = 0.d0
!
! - Coordinates of nodes for SB7 elements: I, J, K, L, M, N
!
    Xi = geom(1)
    Yi = geom(2)
    Zi = geom(3)
    Xj = geom(4)
    Yj = geom(5)
    Zj = geom(6)
    Xk = geom(7)
    Yk = geom(8)
    Zk = geom(9)
    Xl = geom(10)
    Yl = geom(11)
    Zl = geom(12)
    Xm = geom(13)
    Ym = geom(14)
    Zm = geom(15)
    Xn = geom(16)
    Yn = geom(17)
    Zn = geom(18)
!
! - Face for THETA_4
!
    call sshGetFaceSB7(Xi, Yi, Zi, Xj, Yj, Zj,&
                       Xm, Ym, Zm, Xl, Yl, Zl,&
                       xS, yS, zS,&
                       xQ, yQ, zQ,&
                       xR, yR, zR,&
                       si, qi, sj, qj,&
                       sl, ql, sm, qm,&
                       faceArea)
!
! - Compute THETA_4
!
    Pf=(sl-sj)/faceArea/2.d0
    T(1,1)=Pf*xR
    T(1,2)=Pf*yR
    T(1,3)=Pf*zR
    T(1,13)=-Pf*xR
    T(1,14)=-Pf*yR
    T(1,15)=-Pf*zR
    Qf=(si-sm)/faceArea/2.d0
    T(1,4)=Qf*xR
    T(1,5)=Qf*yR
    T(1,6)=Qf*zR
    T(1,10)=-Qf*xR
    T(1,11)=-Qf*yR
    T(1,12)=-Qf*zR
!
! - Face for THETA_5
!
    call sshGetFaceSB7(Xj, Yj, Zj, Xk, Yk, Zk,&
                       Xn, Yn, Zn, Xm, Ym, Zm,&
                       xS, yS, zS,&
                       xQ, yQ, zQ,&
                       xR, yR, zR,&
                       sj, qj, sk, qk,&
                       sm, qm, sn, qn,&
                       faceArea)
!
! - Compute THETA_5
!
    Pf=(sm-sk)/faceArea/2.d0
    T(2,4)=Pf*xR
    T(2,5)=Pf*yR
    T(2,6)=Pf*zR
    T(2,16)=-Pf*xR
    T(2,17)=-Pf*yR
    T(2,18)=-Pf*zR
    Qf=(sj-sn)/faceArea/2.d0
    T(2,7)=Qf*xR
    T(2,8)=Qf*yR
    T(2,9)=Qf*zR
    T(2,13)=-Qf*xR
    T(2,14)=-Qf*yR
    T(2,15)=-Qf*zR
!
! - Face for THETA_6
!
    call sshGetFaceSB7(Xk, Yk, Zk, Xi, Yi, Zi,&
                       Xl, Yl, Zl, Xn, Yn, Zn,&
                       xS, yS, zS,&
                       xQ, yQ, zQ,&
                       xR, yR, zR,&
                       sk, qk, si, qi,&
                       sn, qn, sl, ql,&
                       faceArea)
!
! - Compute THETA_6
!
    Pf=(sn-si)/faceArea/2.d0
    T(3,7)=Pf*xR
    T(3,8)=Pf*yR
    T(3,9)=Pf*zR
    T(3,10)=-Pf*xR
    T(3,11)=-Pf*yR
    T(3,12)=-Pf*zR
    Qf=(sk-sl)/faceArea/2.d0
    T(3,1)=Qf*xR
    T(3,2)=Qf*yR
    T(3,3)=Qf*zR
    T(3,16)=-Qf*xR
    T(3,17)=-Qf*yR
    T(3,18)=-Qf*zR
!
end subroutine
