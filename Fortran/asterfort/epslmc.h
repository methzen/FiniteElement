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
    subroutine epslmc(nno  , ndim  , nbsig,&
                      npg  , ipoids, ivf  ,&
                      idfde, xyz   , depl ,&
                      epsi )
        integer     , intent(in)    :: nno
        integer     , intent(in)    :: ndim
        integer     , intent(in)    :: nbsig
        integer     , intent(in)    :: npg
        integer     , intent(in)    :: ipoids
        integer     , intent(in)    :: ivf
        integer     , intent(in)    :: idfde
        real(kind=8), intent(in)    :: xyz(1)
        real(kind=8), intent(in)    :: depl(1)
        real(kind=8), intent(inout) :: epsi(1)
    end subroutine epslmc
end interface
