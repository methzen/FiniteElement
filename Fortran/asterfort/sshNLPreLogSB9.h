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
    subroutine sshNLPreLogSB9(lVect   , lgpg    , vim ,&
                              epsgPrev, epsgCurr, &
                              epslPrev, epslIncr,&
                              gn      , lamb    , logl,&
                              tPrev   , iret)
        aster_logical, intent(in) :: lVect
        integer, intent(in) :: lgpg
        real(kind=8), intent(in) :: vim(lgpg)
        real(kind=8), intent(in) :: epsgCurr(6), epsgPrev(6)
        real(kind=8), intent(out) :: epslPrev(6), epslIncr(6)
        real(kind=8), intent(out) :: gn(3, 3), lamb(3), logl(3)
        real(kind=8), intent(out) :: tPrev(6)
        integer, intent(out) :: iret
    end subroutine sshNLPreLogSB9
end interface
