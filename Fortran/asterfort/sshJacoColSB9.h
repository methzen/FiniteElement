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
    subroutine sshJacoColSB9(XI ,&
                             J10, J1ETA, J1ZETA, J1ETAZETA,&
                             J20, J2XI , J2ZETA, J2XIZETA ,&
                             J30, J3ETA, J3XI  , J3XIETA  ,&
                             J1 , J2   , J3    )
        real(kind=8), intent(in) :: XI(3)
        real(kind=8), intent(in) :: J10(3), J1ETA(3), J1ZETA(3), J1ETAZETA(3)
        real(kind=8), intent(in) :: J20(3), J2XI(3), J2ZETA(3), J2XIZETA(3)
        real(kind=8), intent(in) :: J30(3), J3ETA(3), J3XI(3), J3XIETA(3)
        real(kind=8), intent(out) :: J1(3), J2(3), J3(3)
    end subroutine sshJacoColSB9
end interface
