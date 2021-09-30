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
subroutine sshNLStabForcSB9(sigm, epsi, dsidep, Ueff)
!
implicit none
!
real(kind=8), intent(in) :: sigm(6), epsi(6),dsidep(6,6)
real(kind=8), intent(out) :: Ueff
!
! --------------------------------------------------------------------------------------------------
!
! Solid-shell element - SB9
!
! Compute effective shear modulus for stabilization
!
! --------------------------------------------------------------------------------------------------
!
! In  sigm             : Cauchy stresses
! In  epsi             : Green-Lagrange strains
! In  dsidep           : jacobian matrix of behaviour (dSigm/dEpsi)
! Out Ueff             : effective shear modulus
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: alpha1, alpha2, Sdev(6), Edev(6)
    real(kind=8) :: aux1, aux2
    integer :: i
    real(kind=8), parameter :: IDEN(6) = (/1.d0, 1.d0, 1.d0, 0.d0, 0.d0, 0.d0/)
!
! --------------------------------------------------------------------------------------------------
!
    Ueff = 0.d0
!
! - Volumic part
!
    alpha1 = (sigm(1)+sigm(2)+sigm(3))/3.d0
    alpha2 = (epsi(1)+epsi(2)+epsi(3))/3.d0
!
! - Deviatoric part
!
    do i = 1, 6
        Sdev(i) = sigm(i) - alpha1*IDEN(i)
        Edev(i) = epsi(i) - alpha2*IDEN(i)
    end do
!
    aux1 = Sdev(1)*Sdev(1)+Sdev(2)*Sdev(2)+Sdev(3)*Sdev(3)+&
           2.d0*(Sdev(4)*Sdev(4)+Sdev(5)*Sdev(5)+Sdev(6)*Sdev(6))
!
    aux2 = Edev(1)*Edev(1)+Edev(2)*Edev(2)+Edev(3)*Edev(3)+&
           (Edev(4)*Edev(4)+Edev(5)*Edev(5)+Edev(6)*Edev(6))/2.d0

    if (abs(aux2).le.0.001d0) then
        Ueff=(dsidep(4,4)+dsidep(5,5)+dsidep(6,6))/3.d0
    else
        Ueff = 0.5d0*sqrt(aux1/aux2)
    endif
!
end subroutine
