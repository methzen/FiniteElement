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
subroutine sshGradMatrPinchSB7(R, bz, Bp)
!
implicit none
!
real(kind=8), intent(in) :: R(3,3), bz(6)
real(kind=8), intent(out) :: Bp(18)
!
! --------------------------------------------------------------------------------------------------
!
! Solid-shell element - SB7
!
! Compute pinching terms of the gradient B matrix
!
! --------------------------------------------------------------------------------------------------
!
! In  R                : rotation matrix for frame of mid_edge triangle
! In  bz               : gradient matrix in local covariant basis - Z coordinate
! Out Bp               : pinch part of the gradient B matrix
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: triaTau1X, triaTau1Y, triaTau1Z
!
! --------------------------------------------------------------------------------------------------
!
    Bp = 0.d0
!
! - Get basis for mid-edge triangle
!
    triaTau1X = R(3,1)
    triaTau1Y = R(3,2)
    triaTau1Z = R(3,3)
!
! - Compute pinch part of the gradient B matrix
!
    Bp(1 )=bz(1)*triaTau1X
    Bp(2 )=bz(1)*triaTau1Y
    Bp(3 )=bz(1)*triaTau1Z
    Bp(4 )=bz(2)*triaTau1X
    Bp(5 )=bz(2)*triaTau1Y
    Bp(6 )=bz(2)*triaTau1Z
    Bp(7 )=bz(3)*triaTau1X
    Bp(8 )=bz(3)*triaTau1Y
    Bp(9 )=bz(3)*triaTau1Z
    Bp(10)=bz(4)*triaTau1X
    Bp(11)=bz(4)*triaTau1Y
    Bp(12)=bz(4)*triaTau1Z
    Bp(13)=bz(5)*triaTau1X
    Bp(14)=bz(5)*triaTau1Y
    Bp(15)=bz(5)*triaTau1Z
    Bp(16)=bz(6)*triaTau1X
    Bp(17)=bz(6)*triaTau1Y
    Bp(18)=bz(6)*triaTau1Z
!
end subroutine
