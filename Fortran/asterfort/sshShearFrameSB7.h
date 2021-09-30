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
    subroutine sshShearFrameSB7(XiL , YiL , ZiL , XjL, YjL, ZjL,&
                                XkL , YkL , ZkL , XlL, YlL, ZlL,&
                                XmL , YmL , ZmL , XnL, YnL, ZnL,&
                                g304, g305, g306,&
                                V1  , V2  , V3  ,&
                                h1  , h2  , h3)
        real(kind=8), intent(in) :: XiL, YiL, ZiL, XjL, YjL, ZjL, XkL, YkL, ZkL
        real(kind=8), intent(in) :: XlL, YlL, ZlL, XmL, YmL, ZmL, XnL, YnL, ZnL
        real(kind=8), intent(out) :: g304(3), g305(3), g306(3)
        real(kind=8), intent(out) :: V1(3,3), V2(3,3), V3(3,3)
        real(kind=8), intent(out) :: h1, h2, h3
    end subroutine sshShearFrameSB7
end interface
