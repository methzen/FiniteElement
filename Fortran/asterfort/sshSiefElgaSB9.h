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
    subroutine sshSiefElgaSB9(npg    , nbsig   ,&
                              nb_node, nb_dof  ,&
                              j_mater, jv_coopg,&
                              geom   , disp    ,&
                              sigm)
        integer, intent(in) :: npg, nb_node, nb_dof, nbsig
        integer, intent(in) :: jv_coopg, j_mater
        real(kind=8), intent(in) :: geom(3*nb_node), disp(nb_dof)
        real(kind=8), intent(out) :: sigm(nbsig*npg)
    end subroutine sshSiefElgaSB9
end interface
