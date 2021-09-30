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
    subroutine cm07ma(mesh_in     , mesh_out    ,&
                      nb_node_mesh, nb_elem_mesh, nb_penta6 ,&
                      nb_elem     , list_elem   ,&
                      v_type_in   , v_type_out)
        character(len=8), intent(in) :: mesh_out, mesh_in
        integer, intent(in) :: nb_node_mesh, nb_elem_mesh, nb_penta6
        integer, intent(in) :: nb_elem, list_elem(nb_elem)
        integer, pointer :: v_type_in(:), v_type_out(:)
    end subroutine cm07ma
end interface
