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
subroutine sshDispDeriCovaSB9(nb_node, nb_dof, disp  ,&
                              D10    , D1ETA , D1ZETA, D1ETAZETA,&
                              D20    , D2XI  , D2ZETA, D2XIZETA ,&
                              D30    , D3ETA , D3XI  , D3XIETA)
!
implicit none
!
#include "asterfort/assert.h"
!
integer, intent(in) :: nb_node, nb_dof
real(kind=8), intent(in) :: disp(nb_dof)
real(kind=8), intent(out) :: D10(3), D1ETA(3), D1ZETA(3), D1ETAZETA(3)
real(kind=8), intent(out) :: D20(3), D2XI(3), D2ZETA(3), D2XIZETA(3)
real(kind=8), intent(out) :: D30(3), D3ETA(3), D3XI(3), D3XIETA(3)
!
! --------------------------------------------------------------------------------------------------
!
! Solid-shell element - SB9
!
! Compute the decomposed derivatives of the displacements in the covariant base
!
! --------------------------------------------------------------------------------------------------
!
! Decomposition of strains:
!   Dc  = Dc0 + xi.DcXI + eta.DcETA + zeta.DcZETA +
!         xi.eta.DcXIETA + eta.zeta.DcETAZETA + xi.zeta.DcXIZETA
!
!   Dc1 = (g1 +        + h1.eta  + h3.zeta + h4.eta.zeta) . U
!   Dc2 = (g2 + h1.xi  +           h2.zeta + h4.xi.zeta ) . U
!   Dc3 = (g3 + h3.xi  + h2.eta            + h4.xi.eta  ) . U
!
! With:
!   Dci = (d X_j / d xi_i) e_j = dN/dxi_i . X <=> J = [J1 J2 J3]
!
! In  nb_node          : number of nodes of element without pinch node(s)
! In  nb_dof           : number of dof
! In  disp             : displacements of element (U)
! Out D10              : decomposed derivatives of the disp. in the covariant base - g1 . U
! Out D20              : decomposed derivatives of the disp. in the covariant base - g2 . U
! Out D30              : decomposed derivatives of the disp. in the covariant base - g3 . U
! Out D1ETA            : decomposed derivatives of the disp. in the covariant base - h1 . U
! Out D1ZETA           : decomposed derivatives of the disp. in the covariant base - h3 . U
! Out D1ETAZETA        : decomposed derivatives of the disp. in the covariant base - h4 . U
! Out D2XI             : decomposed derivatives of the disp. in the covariant base - h1 . U
! Out D2ZETA           : decomposed derivatives of the disp. in the covariant base - h2 . U
! Out D2XIZETA         : decomposed derivatives of the disp. in the covariant base - h4 . U
! Out D3XI             : decomposed derivatives of the disp. in the covariant base - h3 . U
! Out D3ETA            : decomposed derivatives of the disp. in the covariant base - h2 . U
! Out D3XIETA          : decomposed derivatives of the disp. in the covariant base - h4 . U
!
! --------------------------------------------------------------------------------------------------
!
    integer :: i
    real(kind=8):: U(8), V(8), W(8)
    real(kind=8):: g1(8), g2(8), g3(8)
    real(kind=8):: h1(8), h2(8), h3(8), h4(8)
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(nb_node .eq. 8)
!
    do i = 1, nb_node
        U(i) = disp(3*(i-1)+1)
        V(i) = disp(3*(i-1)+2)
        W(i) = disp(3*(i-1)+3)
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
    D10(1)       = sum(g1*U)
    D10(2)       = sum(g1*V)
    D10(3)       = sum(g1*W)
    D1ETA(1)     = sum(h1*U)
    D1ETA(2)     = sum(h1*V)
    D1ETA(3)     = sum(h1*W)
    D1ZETA(1)    = sum(h3*U)
    D1ZETA(2)    = sum(h3*V)
    D1ZETA(3)    = sum(h3*W)
    D1ETAZETA(1) = sum(h4*U)
    D1ETAZETA(2) = sum(h4*V)
    D1ETAZETA(3) = sum(h4*W)
!
    D20(1)       = sum(g2*U)
    D20(2)       = sum(g2*V)
    D20(3)       = sum(g2*W)
    D2XI(1)      = sum(h1*U)
    D2XI(2)      = sum(h1*V)
    D2XI(3)      = sum(h1*W)
    D2ZETA(1)    = sum(h2*U)
    D2ZETA(2)    = sum(h2*V)
    D2ZETA(3)    = sum(h2*W)
    D2XIZETA(1)  = sum(h4*U)
    D2XIZETA(2)  = sum(h4*V)
    D2XIZETA(3)  = sum(h4*W)
!
    D30(1)       = sum(g3*U)
    D30(2)       = sum(g3*V)
    D30(3)       = sum(g3*W)
    D3ETA(1)     = sum(h2*U)
    D3ETA(2)     = sum(h2*V)
    D3ETA(3)     = sum(h2*W)
    D3XI(1)      = sum(h3*U)
    D3XI(2)      = sum(h3*V)
    D3XI(3)      = sum(h3*W)
    D3XIETA(1)   = sum(h4*U)
    D3XIETA(2)   = sum(h4*V)
    D3XIETA(3)   = sum(h4*W)
!
end subroutine
