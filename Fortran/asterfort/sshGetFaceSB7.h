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
interface
    subroutine sshGetFaceSB7(Xi, Yi, Zi, Xj, Yj, Zj,&
                             Xm, Ym, Zm, Xl, Yl, Zl,&
                             xS, yS, zS,&
                             xQ, yQ, zQ,&
                             xR, yR, zR,&
                             si, qi, sj, qj,&
                             sl, ql, sm, qm,&
                             faceArea)
        real(kind=8), intent(in) :: Xi, Yi, Zi, Xj, Yj, Zj, Xl, Yl, Zl, Xm, Ym, Zm
        real(kind=8), intent(out) :: xR, yR, zR
        real(kind=8), intent(out) :: xS, yS, zS
        real(kind=8), intent(out) :: xQ, yQ, zQ
        real(kind=8), intent(out) :: si, qi, sj, qj, sl, ql, sm, qm
        real(kind=8), intent(out) :: faceArea
    end subroutine sshGetFaceSB7
end interface
