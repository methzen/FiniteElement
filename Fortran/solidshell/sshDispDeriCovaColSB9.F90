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
subroutine sshDispDeriCovaColSB9(XI ,&
                                 D10, D1ETA, D1ZETA, D1ETAZETA,&
                                 D20, D2XI , D2ZETA, D2XIZETA ,&
                                 D30, D3ETA, D3XI  , D3XIETA  ,&
                                 D1 , D2   , D3    )
!
implicit none
!
real(kind=8), intent(in) :: XI(3)
real(kind=8), intent(in) :: D10(3), D1ETA(3), D1ZETA(3), D1ETAZETA(3)
real(kind=8), intent(in) :: D20(3), D2XI(3), D2ZETA(3), D2XIZETA(3)
real(kind=8), intent(in) :: D30(3), D3ETA(3), D3XI(3), D3XIETA(3)
real(kind=8), intent(out) :: D1(3), D2(3), D3(3)
!
! --------------------------------------------------------------------------------------------------
!
! Solid-shell element - SB9
!
! Get columns of the decomposed derivatives of the displacements in the covariant base 
! at current point
!
! --------------------------------------------------------------------------------------------------
!
! Decomposition of gradient:
!   Dc  = Dc0 + xi.DcXI + eta.DcETA + zeta.DcZETA +
!         xi.eta.DcXIETA + eta.zeta.DcETAZETA + xi.zeta.DcXIZETA
!
!   Dc1 = (g1 +        + h1.eta  + h3.zeta + h4.eta.zeta) . U
!   Dc2 = (g2 + h1.xi  +           h2.zeta + h4.xi.zeta ) . U
!   Dc3 = (g3 + h3.xi  + h2.eta            + h4.xi.eta  ) . U
!
! With:
!   Dc = [D1 D2 D3]
!
! In  XI               : parametric coordinates (X)
! In  D10              : decomposed derivatives of the disp. in the covariant base - g1 . U
! In  D20              : decomposed derivatives of the disp. in the covariant base - g2 . U
! In  D30              : decomposed derivatives of the disp. in the covariant base - g3 . U
! In  D1ETA            : decomposed derivatives of the disp. in the covariant base - h1 . U
! In  D1ZETA           : decomposed derivatives of the disp. in the covariant base - h3 . U
! In  D1ETAZETA        : decomposed derivatives of the disp. in the covariant base - h4 . U
! In  D2XI             : decomposed derivatives of the disp. in the covariant base - h1 . U
! In  D2ZETA           : decomposed derivatives of the disp. in the covariant base - h2 . U
! In  D2XIZETA         : decomposed derivatives of the disp. in the covariant base - h4 . U
! In  D3XI             : decomposed derivatives of the disp. in the covariant base - h3 . U
! In  D3ETA            : decomposed derivatives of the disp. in the covariant base - h2 . U
! In  D3XIETA          : decomposed derivatives of the disp. in the covariant base - h4 . U
! Out D1               : first column of the gradient in the covariant base
! Out D2               : second column of the gradient in the covariant base
! Out D3               : third column of the gradient in the covariant base
!
! --------------------------------------------------------------------------------------------------
!
    D1 = D10+XI(2)*D1ETA+XI(3)*D1ZETA+XI(2)*XI(3)*D1ETAZETA
    D2 = D20+XI(1)*D2XI+XI(3)*D2ZETA+XI(1)*XI(3)*D2XIZETA
    D3 = D30+XI(2)*D3ETA+XI(1)*D3XI+XI(1)*XI(2)*D3XIETA
!
end subroutine
