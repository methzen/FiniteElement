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
    subroutine sshDispDeriCovaSB9(nb_node, nb_dof, disp  ,&
                                  D10    , D1ETA , D1ZETA, D1ETAZETA,&
                                  D20    , D2XI  , D2ZETA, D2XIZETA ,&
                                  D30    , D3ETA , D3XI  , D3XIETA)
        integer, intent(in) :: nb_node, nb_dof
        real(kind=8), intent(in) :: disp(nb_dof)
        real(kind=8), intent(out) :: D10(3), D1ETA(3), D1ZETA(3), D1ETAZETA(3)
        real(kind=8), intent(out) :: D20(3), D2XI(3), D2ZETA(3), D2XIZETA(3)
        real(kind=8), intent(out) :: D30(3), D3ETA(3), D3XI(3), D3XIETA(3)
    end subroutine sshDispDeriCovaSB9
end interface
