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
    subroutine sshNLPosLogSB9(lVect , lMatr  , lgpg,&
                              tPrev , tCurr  ,&
                              gn    , lamb   , logl,&
                              dtde  , vip    , &
                              dsidep, pk2Curr)
        aster_logical, intent(in) :: lVect, lMatr
        integer, intent(in) :: lgpg
        real(kind=8), intent(in) :: tPrev(6), tCurr(6)
        real(kind=8), intent(in) :: gn(3, 3), lamb(3), logl(3)
        real(kind=8), intent(in) :: dtde(6, 6)
        real(kind=8), intent(inout) :: vip(lgpg)
        real(kind=8), intent(out) :: dsidep(6, 6), pk2Curr(6)
    end subroutine sshNLPosLogSB9
end interface
