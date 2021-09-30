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
    subroutine sshRigiMatrSB7(npg        , nb_node , nb_dof  ,&
                              jv_geom    , jv_coopg, jv_poids,&
                              hookeMatrix, matrRigi)
        integer, intent(in) :: npg, nb_node, nb_dof
        integer, intent(in) :: jv_geom, jv_poids, jv_coopg
        real(kind=8), intent(in) :: hookeMatrix(6, 6)
        real(kind=8), intent(out) :: matrRigi(25, 25)
    end subroutine sshRigiMatrSB7
end interface
