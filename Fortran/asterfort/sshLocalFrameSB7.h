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
    subroutine sshLocalFrameSB7(geom   ,&
                                R      , geomLocal, h   ,&
                                triaN1x_ , triaN1y_ ,&
                                triaN2x_ , triaN2y_ ,&
                                triaN3x_ , triaN3y_ ,&
                                triaArea_)
        real(kind=8), intent(in)  :: geom(18)
        real(kind=8), intent(out) :: geomLocal(18), R(3,3), h
        real(kind=8), optional, intent(out) :: triaN1x_, triaN1y_
        real(kind=8), optional, intent(out) :: triaN2x_, triaN2y_
        real(kind=8), optional, intent(out) :: triaN3x_, triaN3y_
        real(kind=8), optional, intent(out) :: triaArea_
    end subroutine sshLocalFrameSB7
end interface
