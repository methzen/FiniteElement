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
subroutine sshGradMatrMembSB7(R, bx, by, Bm)
!
implicit none
!
real(kind=8), intent(in) :: R(3,3)
real(kind=8), intent(in) :: bx(6), by(6)
real(kind=8), intent(out) :: Bm(3,18)
!
! --------------------------------------------------------------------------------------------------
!
! Solid-shell element - SB7
!
! Compute membrane terms of the gradient B matrix
!
! --------------------------------------------------------------------------------------------------
!
! In  R                : rotation matrix for frame of mid_edge triangle
! In  bx               : gradient matrix in local covariant basis - X coordinate
! In  by               : gradient matrix in local covariant basis - Y coordinate
! Out Bm               : membrane part of the gradient B matrix
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: triaNormX, triaNormY, triaNormZ
    real(kind=8) :: triaTau2X, triaTau2Y, triaTau2Z
!
! --------------------------------------------------------------------------------------------------
!
    Bm = 0.d0
!
! - Get basis for mid-edge triangle
!
    triaNormX = R(1,1)
    triaNormY = R(1,2)
    triaNormZ = R(1,3)
    triaTau2X = R(2,1)
    triaTau2Y = R(2,2)
    triaTau2Z = R(2,3)
!
! - Compute membrane part of the gradient B matrix
!
    Bm(1,1 )=bx(1)*triaNormX
    Bm(2,1 )=by(1)*triaTau2X
    Bm(3,1 )=bx(1)*triaTau2X+by(1)*triaNormX
    Bm(1,2 )=bx(1)*triaNormY
    Bm(2,2 )=by(1)*triaTau2Y
    Bm(3,2 )=bx(1)*triaTau2Y+by(1)*triaNormY
    Bm(1,3 )=bx(1)*triaNormZ
    Bm(2,3 )=by(1)*triaTau2Z
    Bm(3,3 )=bx(1)*triaTau2Z+by(1)*triaNormZ
    Bm(1,4 )=bx(2)*triaNormX
    Bm(2,4 )=by(2)*triaTau2X
    Bm(3,4 )=bx(2)*triaTau2X+by(2)*triaNormX
    Bm(1,5 )=bx(2)*triaNormY
    Bm(2,5 )=by(2)*triaTau2Y
    Bm(3,5 )=bx(2)*triaTau2Y+by(2)*triaNormY
    Bm(1,6 )=bx(2)*triaNormZ
    Bm(2,6 )=by(2)*triaTau2Z
    Bm(3,6 )=bx(2)*triaTau2Z+by(2)*triaNormZ
    Bm(1,7 )=bx(3)*triaNormX
    Bm(2,7 )=by(3)*triaTau2X
    Bm(3,7 )=bx(3)*triaTau2X+by(3)*triaNormX
    Bm(1,8 )=bx(3)*triaNormY
    Bm(2,8 )=by(3)*triaTau2Y
    Bm(3,8 )=bx(3)*triaTau2Y+by(3)*triaNormY
    Bm(1,9 )=bx(3)*triaNormZ
    Bm(2,9 )=by(3)*triaTau2Z
    Bm(3,9 )=bx(3)*triaTau2Z+by(3)*triaNormZ
    Bm(1,10)=bx(4)*triaNormX
    Bm(2,10)=by(4)*triaTau2X
    Bm(3,10)=bx(4)*triaTau2X+by(4)*triaNormX
    Bm(1,11)=bx(4)*triaNormY
    Bm(2,11)=by(4)*triaTau2Y
    Bm(3,11)=bx(4)*triaTau2Y+by(4)*triaNormY
    Bm(1,12)=bx(4)*triaNormZ
    Bm(2,12)=by(4)*triaTau2Z
    Bm(3,12)=bx(4)*triaTau2Z+by(4)*triaNormZ
    Bm(1,13)=bx(5)*triaNormX
    Bm(2,13)=by(5)*triaTau2X
    Bm(3,13)=bx(5)*triaTau2X+by(5)*triaNormX
    Bm(1,14)=bx(5)*triaNormY
    Bm(2,14)=by(5)*triaTau2Y
    Bm(3,14)=bx(5)*triaTau2Y+by(5)*triaNormY
    Bm(1,15)=bx(5)*triaNormZ
    Bm(2,15)=by(5)*triaTau2Z
    Bm(3,15)=bx(5)*triaTau2Z+by(5)*triaNormZ
    Bm(1,16)=bx(6)*triaNormX
    Bm(2,16)=by(6)*triaTau2X
    Bm(3,16)=bx(6)*triaTau2X+by(6)*triaNormX
    Bm(1,17)=bx(6)*triaNormY
    Bm(2,17)=by(6)*triaTau2Y
    Bm(3,17)=bx(6)*triaTau2Y+by(6)*triaNormY
    Bm(1,18)=bx(6)*triaNormZ
    Bm(2,18)=by(6)*triaTau2Z
    Bm(3,18)=bx(6)*triaTau2Z+by(6)*triaNormZ
!
end subroutine
