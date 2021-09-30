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
subroutine sshTMatrDecoSB9(nb_node, geom, TXI_, TETA_, TZETA_)
!
implicit none
!
#include "asterfort/assert.h"
#include "asterfort/matinv.h"
#include "asterfort/sshJaco.h"
#include "asterfort/sshTMatrCmpDecoSB9.h"
!
integer, intent(in) :: nb_node
real(kind=8), intent(in) :: geom(3*nb_node)
real(kind=8), optional, intent(out) :: TXI_(6,6), TETA_(6,6), TZETA_(6,6)
!
! --------------------------------------------------------------------------------------------------
!
! Solid-shell element - SB9
!
! Compute the decomposition of the T matrix (matrix relating the covariant and cartesian frames)
!
! --------------------------------------------------------------------------------------------------
!
! T = [TXI TETA TZETA] = T0 + xi.TXI + eta.TETA + zeta.TZETA
!
! In  nb_node          : number of nodes of element without pinch node(s)
! In  geom             : coordinates of element
! Out TXI              : xi components for T matrix
! Out TETA             : eta components for T matrix
! Out TZETA            : zeta components for T matrix
!
! --------------------------------------------------------------------------------------------------
!
    integer :: i
    real(kind=8) :: X(8), Y(8), Z(8), h1(8), h2(8), h3(8)
    real(kind=8) :: J0(3,3), JXI(3,3), JETA(3,3), JZETA(3,3)
    real(kind=8) :: JXIJ0(3,3), JETAJ0(3,3)
    real(kind=8) :: JZETAJ0(3,3), A(3,3), B(3,3), C(3,3)
    real(kind=8) :: Jac(3,3), det
    real(kind=8) :: TXI(6,6), TETA(6,6), TZETA(6,6)
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(nb_node .eq. 8)
    TXI   = 0.d0
    TETA  = 0.d0
    TZETA = 0.d0
!
    do i = 1, 8
        X(i) = geom(3*(i-1)+1)
        Y(i) = geom(3*(i-1)+2)
        Z(i) = geom(3*(i-1)+3)
    enddo
!
    h1=(/1.d0, -1.d0,  1.d0, -1.d0,  1.d0, -1.d0, 1.d0, -1.d0/)/8.d0
    h2=(/1.d0,  1.d0, -1.d0, -1.d0, -1.d0, -1.d0, 1.d0,  1.d0/)/8.d0
    h3=(/1.d0, -1.d0, -1.d0,  1.d0, -1.d0,  1.d0, 1.d0, -1.d0/)/8.d0
!
! - Compute jacobian matrix at center of element
!
    call sshJaco(nb_node, geom, (/0.d0, 0.d0, 0.d0 /), Jac)
!
! - Compute inverse jacobian matrix
!
    call matinv('S', 3, Jac, J0, det)
!
    A        = 0.d0
    B        = 0.d0
    C        = 0.d0
    JXIJ0    = 0.d0
    JETAJ0   = 0.d0
    JZETAJ0  = 0.d0
!
    JXI(1,1) = 0.d0
    JXI(1,2) = sum(X*h1)
    JXI(1,3) = sum(X*h3)
    JXI(2,1) = 0.d0
    JXI(2,2) = sum(Y*h1)
    JXI(2,3) = sum(Y*h3)
    JXI(3,1) = 0.d0
    JXI(3,2) = sum(Z*h1)
    JXI(3,3) = sum(Z*h3)
!
    JETA(1,1) = sum(X*h1)
    JETA(1,2) = 0.d0
    JETA(1,3) = sum(X*h2)
    JETA(2,1) = sum(Y*h1)
    JETA(2,2) = 0.d0
    JETA(2,3) = sum(Y*h2)
    JETA(3,1) = sum(Z*h1)
    JETA(3,2) = 0.d0
    JETA(3,3) = sum(Z*h2)
!
    JZETA(1,1) = sum(X*h3)
    JZETA(1,2) = sum(X*h2)
    JZETA(1,3) = 0.d0
    JZETA(2,1) = sum(Y*h3)
    JZETA(2,2) = sum(Y*h2)
    JZETA(2,3) = 0.d0
    JZETA(3,1) = sum(Z*h3)
    JZETA(3,2) = sum(Z*h2)
    JZETA(3,3) = 0.d0
!
    JXIJ0   = matmul(JXI,J0)
    JETAJ0  = matmul(JETA,J0)
    JZETAJ0 = matmul(JZETA,J0)
!
    A       = -matmul(J0,JXIJ0)
    B       = -matmul(J0,JETAJ0)
    C       = -matmul(J0,JZETAJ0)
!
    call sshTMatrCmpDecoSB9(J0, A, TXI)
    call sshTMatrCmpDecoSB9(J0, B, TETA)
    call sshTMatrCmpDecoSB9(J0, C, TZETA)
!
    if (present(TXI_)) then
        TXI_ = TXI
    endif
    if (present(TETA_)) then
        TETA_ = TETA
    endif
    if (present(TZETA_)) then
        TZETA_ = TZETA
    endif
!
end subroutine
