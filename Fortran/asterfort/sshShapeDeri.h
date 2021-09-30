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
    subroutine sshShapeDeri(nb_node, XI, dN_dXsi, dN_dEta, dN_dZeta)
        integer, intent(in) :: nb_node
        real(kind=8), intent(in)  :: XI(3)
        real(kind=8), intent(out) :: dN_dXsi(nb_node), dN_dEta(nb_node), dN_dZeta(nb_node)
    end subroutine sshShapeDeri
end interface
