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
    subroutine sshGreenCartSB9(nb_node , nb_dof,&
                               geomInit, disp  ,&
                               EH1     , EH2   ,&
                               EH23    , EH13)
        integer, intent(in) :: nb_node, nb_dof
        real(kind=8), intent(in) :: geomInit(3*nb_node), disp(nb_dof)
        real(kind=8), intent(out) :: EH1(6), EH2(6), EH23(6), EH13(6)
    end subroutine sshGreenCartSB9
end interface
