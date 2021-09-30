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
    subroutine sshGradMatrCovaSB9(nb_node , geom     , Bc0       ,&
                                  BcXI    , BcETA    , BcZETA    ,&
                                  BcXIZETA, BcETAZETA, BcZETAZETA)
        integer, intent(in) :: nb_node
        real(kind=8), intent(in)  :: geom(3*nb_node)
        real(kind=8), intent(out) :: Bc0(6,3*nb_node), BcZETA(6,3*nb_node)
        real(kind=8), intent(out) :: BcZETAZETA(6,3*nb_node), BcXI(6,3*nb_node)
        real(kind=8), intent(out) :: BcETA(6,3*nb_node), BcETAZETA(6,3*nb_node)
        real(kind=8), intent(out) :: BcXIZETA(6,3*nb_node)
    end subroutine sshGradMatrCovaSB9
end interface
