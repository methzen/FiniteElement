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
    subroutine sshNLStabSB7(lMatr   , lVect ,&
                            nb_node , nb_dof,&
                            jacobian, Ueff  ,&
                            dispIncr, BStab ,&
                            matuu   , vectu)
        aster_logical, intent(in) :: lMatr, lVect
        integer, intent(in) :: nb_node, nb_dof
        real(kind=8), intent(in) :: dispIncr(nb_dof)
        real(kind=8), intent(in)  :: Ueff, jacobian
        real(kind=8), intent(in) :: BStab(2, 18)
        real(kind=8), intent(inout) :: vectu(*), matuu(*)
    end subroutine sshNLStabSB7
end interface
