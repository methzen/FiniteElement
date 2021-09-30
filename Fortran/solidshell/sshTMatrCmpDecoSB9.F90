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
! but WITHOUT ANY WARRANTY without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with code_aster.  If not, see <http://www.gnu.org/licenses/>.
! --------------------------------------------------------------------
!
subroutine sshTMatrCmpDecoSB9(J0, A, TXI)
!
implicit none
!
real(kind=8), intent(in)  :: J0(3, 3), A(3, 3)
real(kind=8), intent(out) :: TXI(6, 6)
!
! --------------------------------------------------------------------------------------------------
!
! Solid-shell element - SB9
!
! CALCULE LES JACOBIENS INVERSES 6x6 DECOMPOSES
!
! --------------------------------------------------------------------------------------------------
!
! In  J0: inverse jacobian matrix
!
! --------------------------------------------------------------------------------------------------
!
    TXI(1,1) = 2.d0*J0(1,1)*A(1,1)
    TXI(1,2) = 2.d0*J0(2,1)*A(2,1)
    TXI(1,3) = 2.d0*J0(3,1)*A(3,1)
    TXI(1,4) = J0(1,1)*A(2,1)+J0(2,1)*A(1,1)
    TXI(1,5) = J0(1,1)*A(3,1)+J0(3,1)*A(1,1)
    TXI(1,6) = J0(2,1)*A(3,1)+J0(3,1)*A(2,1)
!
    TXI(2,1) = 2.d0*J0(1,2)*A(1,2)
    TXI(2,2) = 2.d0*J0(2,2)*A(2,2)
    TXI(2,3) = 2.d0*J0(3,2)*A(3,2)
    TXI(2,4) = J0(1,2)*A(2,2)+J0(2,2)*A(1,2)
    TXI(2,5) = J0(1,2)*A(3,2)+J0(3,2)*A(1,2)
    TXI(2,6) = J0(2,2)*A(3,2)+J0(3,2)*A(2,2)
!
    TXI(3,1) = 2.d0*J0(1,3)*A(1,3)
    TXI(3,2) = 2.d0*J0(2,3)*A(2,3)
    TXI(3,3) = 2.d0*J0(3,3)*A(3,3)
    TXI(3,4) = J0(1,3)*A(2,3)+J0(2,3)*A(1,3)
    TXI(3,5) = J0(1,3)*A(3,3)+J0(3,3)*A(1,3)
    TXI(3,6) = J0(2,3)*A(3,3)+J0(3,3)*A(2,3)
!
    TXI(4,1) = 2.d0*(J0(1,1)*A(1,2)+J0(1,2)*A(1,1))
    TXI(4,2) = 2.d0*(J0(2,1)*A(2,2)+J0(2,2)*A(2,1))
    TXI(4,3) = 2.d0*(J0(3,1)*A(3,2)+J0(3,2)*A(3,1))
    TXI(4,4) = (J0(1,2)*A(2,1)+J0(2,1)*A(1,2))+(J0(1,1)*A(2,2)+J0(2,2)*A(1,1))
    TXI(4,5) = (J0(1,2)*A(3,1)+J0(3,1)*A(1,2))+(J0(1,1)*A(3,2)+J0(3,2)*A(1,1))
    TXI(4,6) = (J0(2,2)*A(3,1)+J0(3,1)*A(2,2))+(J0(2,1)*A(3,2)+J0(3,2)*A(2,1))
!
    TXI(5,1) = 2.d0*(J0(1,1)*A(1,3)+J0(1,3)*A(1,1))
    TXI(5,2) = 2.d0*(J0(2,1)*A(2,3)+J0(2,3)*A(2,1))
    TXI(5,3) = 2.d0*(J0(3,1)*A(3,3)+J0(3,3)*A(3,1))
    TXI(5,4) = (J0(1,3)*A(2,1)+J0(2,1)*A(1,3))+(J0(1,1)*A(2,3)+J0(2,3)*A(1,1))
    TXI(5,5) = (J0(1,3)*A(3,1)+J0(3,1)*A(1,3))+(J0(1,1)*A(3,3)+J0(3,3)*A(1,1))
    TXI(5,6) = (J0(2,3)*A(3,1)+J0(3,1)*A(2,3))+(J0(2,1)*A(3,3)+J0(3,3)*A(2,1))
!
    TXI(6,1) = 2.d0*(J0(1,2)*A(1,3)+J0(1,3)*A(1,2))
    TXI(6,2) = 2.d0*(J0(2,2)*A(2,3)+J0(2,3)*A(2,2))
    TXI(6,3) = 2.d0*(J0(3,2)*A(3,3)+J0(3,3)*A(3,2))
    TXI(6,4) = (J0(1,3)*A(2,2)+J0(2,2)*A(1,3))+(J0(1,2)*A(2,3)+J0(2,3)*A(1,2))
    TXI(6,5) = (J0(1,3)*A(3,2)+J0(3,2)*A(1,3))+(J0(1,2)*A(3,3)+J0(3,3)*A(1,2))
    TXI(6,6) = (J0(2,3)*A(3,2)+J0(3,2)*A(2,3))+(J0(2,2)*A(3,3)+J0(3,3)*A(2,2))
!
end subroutine
