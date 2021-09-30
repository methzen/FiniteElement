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
    subroutine sshNLMatr(nbsig   , nb_node, nb_dof,&
                         jacobian, B      , dsidep,&
                         matuu   , kGeom_)
        integer, intent(in) :: nb_node, nb_dof, nbsig
        real(kind=8), intent(in) :: jacobian
        real(kind=8), intent(in) :: dsidep(nbsig, nbsig), B(nbsig, nb_dof)
        real(kind=8), intent(inout) :: matuu(nb_dof*nb_dof)
        real(kind=8), optional, intent(in) :: kGeom_(nb_dof, nb_dof)
    end subroutine sshNLMatr
end interface
