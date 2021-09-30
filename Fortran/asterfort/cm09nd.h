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
    subroutine cm09nd(mesh_in  , mesh_out ,&
                      pref_name, pref_nume, nb_node_mesh,&
                      nb_elem  , list_elem, v_type_elem , v_coor,&
                      nb_hexa8)
        character(len=8), intent(in) :: mesh_in, mesh_out, pref_name
        integer, intent(in) :: pref_nume, nb_node_mesh
        integer, intent(in) :: nb_elem, list_elem(nb_elem)
        integer, pointer :: v_type_elem(:) 
        real(kind=8), intent(inout) :: v_coor(3, *)
        integer, intent(in) :: nb_hexa8
    end subroutine cm09nd
end interface
