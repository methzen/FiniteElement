! --------------------------------------------------------------------
! Copyright (C) 1991 - 2020 - EDF R&J - www.code-aster.org
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
subroutine sshJacoDecoSB9(nb_node, geom ,&
                          J10    , J1ETA, J1ZETA, J1ETAZETA,&
                          J20    , J2XI , J2ZETA, J2XIZETA ,&
                          J30    , J3ETA, J3XI  , J3XIETA)
!
implicit none
!
#include "asterfort/assert.h"
!
integer, intent(in) :: nb_node
real(kind=8), intent(in)  :: geom(3*nb_node)
real(kind=8), intent(out) :: J10(3), J1ETA(3), J1ZETA(3), J1ETAZETA(3)
real(kind=8), intent(out) :: J20(3), J2XI(3), J2ZETA(3), J2XIZETA(3)
real(kind=8), intent(out) :: J30(3), J3ETA(3), J3XI(3), J3XIETA(3)
!
! --------------------------------------------------------------------------------------------------
!
! Solid-shell element - SB9
!
! Compute the decomposed derivatives of the jacobian in the covariant base
!
! --------------------------------------------------------------------------------------------------
!
! Decomposition of jacobian:
!   J  = J0 + xi.JXI + eta.JETA + zeta.JZETA +
!        xi.eta.JXIETA + eta.zeta.JETAZETA + xi.zeta.JXIZETA
!
!   J1 = (g1 +        + h1.eta  + h3.zeta + h4.eta.zeta) . X
!   J2 = (g2 + h1.xi  +           h2.zeta + h4.xi.zeta ) . X
!   J3 = (g3 + h3.xi  + h2.eta            + h4.xi.eta  ) . X
!
! With:
!   Ji = (d X_j / d xi_i) e_j = dN/dxi_i . X <=> J = [J1 J2 J3]
!
! In  nb_node          : number of nodes of element without pinch node(s)
! In  geom             : coordinates of element (X)
! Out J10              : decomposed derivatives of the jacobian in the covariant base - g1 . X
! Out J20              : decomposed derivatives of the jacobian in the covariant base - g2 . X
! Out J30              : decomposed derivatives of the jacobian in the covariant base - g3 . X
! Out J1ETA            : decomposed derivatives of the jacobian in the covariant base - h1 . X
! Out J1ZETA           : decomposed derivatives of the jacobian in the covariant base - h3 . X
! Out J1ETAZETA        : decomposed derivatives of the jacobian in the covariant base - h4 . X
! Out J2XI             : decomposed derivatives of the jacobian in the covariant base - h1 . X
! Out J2ZETA           : decomposed derivatives of the jacobian in the covariant base - h2 . X
! Out J2XIZETA         : decomposed derivatives of the jacobian in the covariant base - h4 . X
! Out J3XI             : decomposed derivatives of the jacobian in the covariant base - h3 . X
! Out J3ETA            : decomposed derivatives of the jacobian in the covariant base - h2 . X
! Out J3XIETA          : decomposed derivatives of the jacobian in the covariant base - h4 . X
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: U(8), V(8), W(8)
    real(kind=8) :: g1(8), g2(8), g3(8)
    real(kind=8) :: h1(8), h2(8), h3(8), h4(8)
    integer :: i
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(nb_node .eq. 8)
    do i = 1, nb_node
        U(i) = geom(3*(i-1)+1)
        V(i) = geom(3*(i-1)+2)
        W(i) = geom(3*(i-1)+3)
    enddo
!
! - Vectors for decomposition
!
    g1=(/-1.d0,  1.d0,  1.d0, -1.d0, -1.d0,  1.d0, 1.d0, -1.d0/)/8.d0
    g2=(/-1.d0, -1.d0,  1.d0,  1.d0, -1.d0, -1.d0, 1.d0,  1.d0/)/8.d0
    g3=(/-1.d0, -1.d0, -1.d0, -1.d0,  1.d0,  1.d0, 1.d0,  1.d0/)/8.d0
    h1=(/ 1.d0, -1.d0,  1.d0, -1.d0,  1.d0, -1.d0, 1.d0, -1.d0/)/8.d0
    h2=(/ 1.d0,  1.d0, -1.d0, -1.d0, -1.d0, -1.d0, 1.d0,  1.d0/)/8.d0
    h3=(/ 1.d0, -1.d0, -1.d0,  1.d0, -1.d0,  1.d0, 1.d0, -1.d0/)/8.d0
    h4=(/-1.d0,  1.d0, -1.d0,  1.d0,  1.d0, -1.d0, 1.d0, -1.d0/)/8.d0
!
    J10(1)       = sum(g1*U)
    J10(2)       = sum(g1*V)
    J10(3)       = sum(g1*W)
    J1ETA(1)     = sum(h1*U)
    J1ETA(2)     = sum(h1*V)
    J1ETA(3)     = sum(h1*W)
    J1ZETA(1)    = sum(h3*U)
    J1ZETA(2)    = sum(h3*V)
    J1ZETA(3)    = sum(h3*W)
    J1ETAZETA(1) = sum(h4*U)
    J1ETAZETA(2) = sum(h4*V)
    J1ETAZETA(3) = sum(h4*W)
!
    J20(1)       = sum(g2*U)
    J20(2)       = sum(g2*V)
    J20(3)       = sum(g2*W)
    J2XI(1)      = sum(h1*U)
    J2XI(2)      = sum(h1*V)
    J2XI(3)      = sum(h1*W)
    J2ZETA(1)    = sum(h2*U)
    J2ZETA(2)    = sum(h2*V)
    J2ZETA(3)    = sum(h2*W)
    J2XIZETA(1)  = sum(h4*U)
    J2XIZETA(2)  = sum(h4*V)
    J2XIZETA(3)  = sum(h4*W)
!
    J30(1)       = sum(g3*U)
    J30(2)       = sum(g3*V)
    J30(3)       = sum(g3*W)
    J3ETA(1)     = sum(h2*U)
    J3ETA(2)     = sum(h2*V)
    J3ETA(3)     = sum(h2*W)
    J3XI(1)      = sum(h3*U)
    J3XI(2)      = sum(h3*V)
    J3XI(3)      = sum(h3*W)
    J3XIETA(1)   = sum(h4*U)
    J3XIETA(2)   = sum(h4*V)
    J3XIETA(3)   = sum(h4*W)
!
end subroutine
