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
subroutine sshGradMatrTwistSB7(geomLocal, R, BStab)
!
implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/sshTwistVectSB7.h"
!
real(kind=8), intent(in) :: geomLocal(18), R(3,3)
real(kind=8), intent(out) :: BStab(2, 18)
!
! --------------------------------------------------------------------------------------------------
!
! Solid-shell element - SB7
!
! Compute stabilization for twist mode
!
! --------------------------------------------------------------------------------------------------
!
! In  geomLocal        : coordinates of element in the frame of mid_edge triangle
! In  R                : rotation matrix for frame of mid_edge triangle
! Out BStab            : pinch part stabilization of the gradient B matrix
!
! --------------------------------------------------------------------------------------------------
!
    integer :: j, k
    real(kind=8) :: Vgamma(6, 2), invJo(3, 3)
!
! --------------------------------------------------------------------------------------------------
!
    BStab = 0.d0
    !WRITE(6,*) 'Stab - GEOMLOCAL: ',sum(geomLocal)
!
! - Compute twist vector for stabilization
!
    call sshTwistVectSB7(geomLocal, invJo, Vgamma)
    !WRITE(6,*) 'Stab - Vgamma: ',sum(Vgamma)
    !WRITE(6,*) 'Stab - invJo: ',sum(invJo)
!
! - Compute 
!
    do j = 1, 6
        k=3*(j-1)
        BStab(1,k+1) = R(3,1)*Vgamma(j,1)
        BStab(1,k+2) = R(3,2)*Vgamma(j,1)
        BStab(1,k+3) = R(3,3)*Vgamma(j,1)
        BStab(2,k+1) = R(3,1)*Vgamma(j,2)
        BStab(2,k+2) = R(3,2)*Vgamma(j,2)
        BStab(2,k+3) = R(3,3)*Vgamma(j,2)
    enddo
    BStab = BStab*invJo(3,3)
!
end subroutine
