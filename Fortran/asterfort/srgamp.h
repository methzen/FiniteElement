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
    subroutine srgamp(val, varv, im, sm, ucrip,&
                      seuilp, vinm, nvi, nbmat, mater, de,&
                      deps, depsv, dgamv, depsp, dgamp, retcom)
        integer :: nbmat
        integer :: val
        integer :: varv
        integer :: nvi
        real(kind=8) :: im
        real(kind=8) :: sm(6)
        real(kind=8) :: ucrip
        real(kind=8) :: seuilp
        real(kind=8) :: vinm(nvi)
        real(kind=8) :: mater(nbmat, 2)
        real(kind=8) :: de(6, 6)
        real(kind=8) :: deps(6)
        real(kind=8) :: depsv(6)
        real(kind=8) :: dgamv
        real(kind=8) :: depsp(6)
        real(kind=8) :: dgamp
        integer :: retcom
    end subroutine srgamp
end interface
