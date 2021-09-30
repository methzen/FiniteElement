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
    subroutine sshBTSigmSB9(npg      , nbsig   ,&
                            nb_node  , nb_dof  ,&
                            jv_coopg , jv_poids,&
                            geom     , sigm    , det      ,&
                            B0       , BZETA   , BZETAZETA,&
                            btSigm)
        integer, intent(in) :: npg, nbsig, nb_node, nb_dof
        integer, intent(in) :: jv_coopg, jv_poids
        real(kind=8), intent(in) :: det
        real(kind=8), intent(in) :: B0(6,24), BZETA(6,24), BZETAZETA(6,24)
        real(kind=8), intent(in) :: geom(3*nb_node), sigm(nbsig*npg)
        real(kind=8), intent(out) :: btSigm(nb_dof)
    end subroutine sshBTSigmSB9
end interface
