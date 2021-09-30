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
    subroutine sshNLSmallSB9(typmod   , option  , fami    ,&
                             npg      , nb_node , nb_dof  , lgpg   ,&
                             jv_poids , jv_coopg, jv_vf   , jv_dfde,&
                             imate     ,&
                             angl_naut, compor  , carcri  ,&
                             timePrev , timeCurr,&
                             geomInit , dispPrev, dispIncr,&
                             sigm     , vim      ,&
                             sigp     , vip      ,&
                             matuu    , vectu    ,&
                             codret)
        integer, intent(in) :: npg, imate, lgpg, nb_node, nb_dof
        integer, intent(in) :: jv_poids, jv_coopg, jv_vf, jv_dfde
        character(len=*), intent(in) :: fami
        character(len=8), intent(in) :: typmod(*)
        character(len=16), intent(in) :: option, compor(*)
        real(kind=8), intent(in) :: carcri(*)
        real(kind=8), intent(in) :: timePrev, timeCurr
        real(kind=8), intent(in) :: geomInit(3*nb_node), dispPrev(nb_dof), dispIncr(nb_dof)
        real(kind=8), intent(in) :: angl_naut(*)
        real(kind=8), intent(in) :: sigm(6,npg), vim(lgpg,npg)
        real(kind=8), intent(out) :: sigp(6,npg), vip(lgpg,npg)
        real(kind=8), intent(out) :: matuu(*), vectu(*)
        integer, intent(out) :: codret
    end subroutine sshNLSmallSB9
end interface
