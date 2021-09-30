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
    subroutine sshGradMatrBendSB7(nb_node, geom, R,&
                                  triaN1x , triaN1y,&
                                  triaN2x , triaN2y,&
                                  triaN3x , triaN3y,&
                                  triaArea, Bb)
        integer, intent(in) :: nb_node
        real(kind=8), intent(in) :: geom(3*nb_node), R(3,3)
        real(kind=8), intent(in) :: triaN1x, triaN1y, triaN2x, triaN2y, triaN3x, triaN3y
        real(kind=8), intent(in) :: triaArea
        real(kind=8), intent(out) :: Bb(3,18)
    end subroutine sshGradMatrBendSB7
end interface
