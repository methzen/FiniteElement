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
subroutine sshJacoColSB9(XI ,&
                         J10, J1ETA, J1ZETA, J1ETAZETA,&
                         J20, J2XI , J2ZETA, J2XIZETA ,&
                         J30, J3ETA, J3XI  , J3XIETA  ,&
                         J1 , J2   , J3    )
!
implicit none
!
real(kind=8), intent(in) :: XI(3)
real(kind=8), intent(in) :: J10(3), J1ETA(3), J1ZETA(3), J1ETAZETA(3)
real(kind=8), intent(in) :: J20(3), J2XI(3), J2ZETA(3), J2XIZETA(3)
real(kind=8), intent(in) :: J30(3), J3ETA(3), J3XI(3), J3XIETA(3)
real(kind=8), intent(out) :: J1(3), J2(3), J3(3)
!
! --------------------------------------------------------------------------------------------------
!
! Solid-shell element - SB9
!
! Get columns of the jacobian in the covariant base at current point
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
! In  XI               : parametric coordinates of current point (X)
! In  J10              : decomposed derivatives of the jacobian in the covariant base - g1 . X
! In  J20              : decomposed derivatives of the jacobian in the covariant base - g2 . X
! In  J30              : decomposed derivatives of the jacobian in the covariant base - g3 . X
! In  J1ETA            : decomposed derivatives of the jacobian in the covariant base - h1 . X
! In  J1ZETA           : decomposed derivatives of the jacobian in the covariant base - h3 . X
! In  J1ETAZETA        : decomposed derivatives of the jacobian in the covariant base - h4 . X
! In  J2XI             : decomposed derivatives of the jacobian in the covariant base - h1 . X
! In  J2ZETA           : decomposed derivatives of the jacobian in the covariant base - h2 . X
! In  J2XIZETA         : decomposed derivatives of the jacobian in the covariant base - h4 . X
! In  J3XI             : decomposed derivatives of the jacobian in the covariant base - h3 . X
! In  J3ETA            : decomposed derivatives of the jacobian in the covariant base - h2 . X
! In  J3XIETA          : decomposed derivatives of the jacobian in the covariant base - h4 . X
! Out J1               : first column of the jacobian in the covariant base
! Out J2               : second column of the jacobian in the covariant base
! Out J3               : third column of the jacobian in the covariant base
!
! --------------------------------------------------------------------------------------------------
!
    J1 = J10+XI(2)*J1ETA+XI(3)*J1ZETA+XI(2)*XI(3)*J1ETAZETA
    J2 = J20+XI(1)*J2XI +XI(3)*J2ZETA+XI(1)*XI(3)*J2XIZETA
    J3 = J30+XI(2)*J3ETA+XI(1)*J3XI  +XI(1)*XI(2)*J3XIETA
!
end subroutine
