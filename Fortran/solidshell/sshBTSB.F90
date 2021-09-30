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
! aslint: disable=W1306
!
subroutine sshBTSB(S, nls, ncb, B, tBSB)
!
implicit none
!
integer, intent(in) :: nls, ncb
real(kind=8), intent(in) :: S(nls,nls), B(nls,ncb)
real(kind=8), intent(out) :: tBSB(ncb,ncb)
!
! --------------------------------------------------------------------------------------------------
!
! Solid-shell element
!
! Compute product tBSB
!
! --------------------------------------------------------------------------------------------------
!
    integer :: i, j, k
    real(kind=8) :: aux, SB(nls,ncb)
!
! --------------------------------------------------------------------------------------------------
!
    tBSB = 0.d0
    do i = 1, nls
        do j = 1, ncb
            aux = 0.0
            do k = 1, nls
                aux = aux+S(i,k)*B(k,j)
                SB(i,j)=aux
            end do
        end do
    end do
!
    do i = 1, ncb
        do j = 1, ncb
            aux = 0.0
            do k = 1, nls
                aux = aux+B(k,i)*SB(k,j)
                tBSB(i,j) = aux
            end do
        end do
    end do
!
end subroutine
