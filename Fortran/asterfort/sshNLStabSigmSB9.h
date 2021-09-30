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
    subroutine sshNLStabSigmSB9(lLarge     ,&
                                nb_node    , nb_dof    ,&
                                geomInit   , disp      , Ueff,&
                                BXIdev     , BETAdev   ,&
                                BETAZETAdev, BXIZETAdev,&
                                SHXI       , SHETA     ,&
                                SHETAZETA  , SHXIZETA)
        aster_logical, intent(in) :: lLarge
        integer, intent(in) :: nb_node, nb_dof
        real(kind=8), intent(in)  :: geomInit(3*nb_node), disp(nb_dof), Ueff
        real(kind=8), intent(in)  :: BXIdev(6,3*nb_node), BETAdev(6,3*nb_node)
        real(kind=8), intent(in)  :: BXIZETAdev(6,3*nb_node), BETAZETAdev(6,3*nb_node)
        real(kind=8), intent(out) :: SHXI(6), SHETA(6), SHETAZETA(6), SHXIZETA(6)
    end subroutine sshNLStabSigmSB9
end interface
