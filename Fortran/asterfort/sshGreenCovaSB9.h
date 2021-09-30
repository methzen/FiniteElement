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
    subroutine sshGreenCovaSB9(nb_node  , nb_dof ,&
                               geomInit , disp   ,&
                               Ec0_     , EcZETA_, EcZETAZETA_,&
                               EcXI_    , EcETA_ , EcETAZETA_,&
                               EcXIZETA_)
        integer, intent(in) :: nb_node, nb_dof
        real(kind=8), intent(in) :: geomInit(3*nb_node), disp(nb_dof)
        real(kind=8), optional, intent(out) :: Ec0_(6), EcZETA_(6), EcZETAZETA_(6)
        real(kind=8), optional, intent(out) :: EcXI_(6), EcETA_(6), EcETAZETA_(6), EcXIZETA_(6)
    end subroutine sshGreenCovaSB9
end interface
