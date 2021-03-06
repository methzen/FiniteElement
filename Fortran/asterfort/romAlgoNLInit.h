! --------------------------------------------------------------------
! Copyright (C) 1991 - 2017 - EDF R&D - www.code-aster.org
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
!
interface
    subroutine romAlgoNLInit(phenom        , model, mesh, nume_dof, result, ds_algorom,&
                             l_line_search_)
        use Rom_Datastructure_type
        character(len=4), intent(in) :: phenom
        character(len=24), intent(in) :: model
        character(len=8), intent(in) :: mesh
        character(len=24), intent(in) :: nume_dof
        character(len=8), intent(in) :: result
        type(ROM_DS_AlgoPara), intent(inout) :: ds_algorom
        aster_logical, intent(in), optional :: l_line_search_
    end subroutine romAlgoNLInit
end interface
