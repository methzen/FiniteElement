! --------------------------------------------------------------------
! Copyright (C) 1991 - 2019 - EDF R&D - www.code-aster.org
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
    subroutine cont_init(mesh  , model    , ds_contact, nume_inst     , ds_measure,&
                         sddyna, hval_incr, sdnume    , list_func_acti)
        use NonLin_Datastructure_type
        character(len=8), intent(in) :: mesh
        character(len=24), intent(in) :: model
        type(NL_DS_Contact), intent(inout) :: ds_contact
        integer, intent(in) :: nume_inst
        type(NL_DS_Measure), intent(inout) :: ds_measure
        character(len=19), intent(in) :: hval_incr(*)
        character(len=19), intent(in) :: sddyna
        integer, intent(in) :: list_func_acti(*)
        character(len=19), intent(in) :: sdnume
    end subroutine cont_init
end interface
