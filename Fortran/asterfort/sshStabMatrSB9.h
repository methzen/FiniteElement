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
    subroutine sshStabMatrSB9(nb_node, Ueff    ,&
                              BXI    , BETA    ,&
                              BXIZETA, BETAZETA,&
                              SXI    , SETA    ,&
                              SXIZETA, SETAZETA)
        integer, intent(in) :: nb_node
        real(kind=8), intent(in) :: Ueff
        real(kind=8), intent(in) :: BXI(6,3*nb_node), BETA(6,3*nb_node)
        real(kind=8), intent(in) :: BXIZETA(6,3*nb_node), BETAZETA(6,3*nb_node)
        real(kind=8), intent(out) :: SXI(3*nb_node,3*nb_node), SETA(3*nb_node,3*nb_node)
        real(kind=8), intent(out) :: SETAZETA(3*nb_node,3*nb_node), SXIZETA(3*nb_node,3*nb_node)
    end subroutine sshStabMatrSB9
end interface
