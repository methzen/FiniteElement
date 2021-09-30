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
subroutine sshGradMatrCovaSB9(nb_node , geom     , Bc0       ,&
                              BcXI    , BcETA    , BcZETA    ,&
                              BcXIZETA, BcETAZETA, BcZETAZETA)
!
implicit none
!
#include "asterfort/assert.h"
#include "asterfort/sshJacoDecoSB9.h"
#include "asterfort/sshJacoColSB9.h"
#include "asterfort/sshShapeDeri.h"
!
integer, intent(in) :: nb_node
real(kind=8), intent(in)  :: geom(3*nb_node)
real(kind=8), intent(out) :: Bc0(6,3*nb_node), BcZETA(6,3*nb_node)
real(kind=8), intent(out) :: BcZETAZETA(6,3*nb_node), BcXI(6,3*nb_node)
real(kind=8), intent(out) :: BcETA(6,3*nb_node), BcETAZETA(6,3*nb_node)
real(kind=8), intent(out) :: BcXIZETA(6,3*nb_node)
!
! --------------------------------------------------------------------------------------------------
!
! Solid-shell element - SB9
!
! Compute gradient matrix in covariant basis
!
! --------------------------------------------------------------------------------------------------
!
! In  nb_node          : number of nodes of element without pinch node(s)
! In  geom             : coordinates of element
! Out Bc0              : constant part of gradient matrix (covariant basis)
! Out BcXI             : xi part of gradient matrix (covariant basis)
! Out BcETA            : eta part of gradient matrix (covariant basis)
! Out BcZETA           : zeta part of gradient matrix (covariant basis)
! Out BcXIZETA         : xi x zeta part of gradient matrix (covariant basis)
! Out BcETAZETA        : eta x zeta part of gradient matrix (covariant basis)
! Out BcZETAZETA       : zeta x zeta part of gradient matrix (covariant basis)
!
! --------------------------------------------------------------------------------------------------
!
    integer :: i, EH, AD, JM, C1, C2
    real(kind=8) :: g1(8), g2(8), g3(8)
    real(kind=8) :: h1(8), h2(8), h3(8), h4(8)
    real(kind=8) :: aux13(3), aux23(3), aux12(3), aux11(3), aux22(3), aux33(3)
    real(kind=8) :: dN_dXsi(8), dN_dEta(8), dN_dZeta(8)
    real(kind=8) :: J10(3), J1ETA(3), J1ZETA(3), J1ETAZETA(3)
    real(kind=8) :: J20(3), J2XI(3), J2ZETA(3), J2XIZETA(3)
    real(kind=8) :: J30(3), J3ETA(3), J3XI(3), J3XIETA(3)
    real(kind=8) :: J1(3), J2(3), J3(3)
    real(kind=8) :: XIEH(4,3), XIAD(4,3), XIJM(4,3)
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(nb_node .eq. 8)
!
    Bc0        = 0.d0
    BcZETA     = 0.d0
    BcZETAZETA = 0.d0
    BcXI       = 0.d0
    BcETA      = 0.d0
    BcETAZETA  = 0.d0
    BcXIZETA   = 0.d0
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
    XIEH(1,1) = -1.d0
    XIEH(1,2) = +0.d0
    XIEH(1,3) = -1.d0
    XIEH(2,1) = +1.d0
    XIEH(2,2) = 0.d0
    XIEH(2,3) = -1.d0
    XIEH(3,1) = +1.d0
    XIEH(3,2) = 0.d0
    XIEH(3,3) = +1.d0
    XIEH(4,1) = -1.d0
    XIEH(4,2) = 0.d0
    XIEH(4,3) = +1.d0
!
    XIAD(1,1) = -1.d0
    XIAD(1,2) = -1.d0
    XIAD(1,3) = 0.d0
    XIAD(2,1) = +1.d0
    XIAD(2,2) = -1.d0
    XIAD(2,3) = 0.d0
    XIAD(3,1) = +1.d0
    XIAD(3,2) = +1.d0
    XIAD(3,3) = 0.d0
    XIAD(4,1) = -1.d0
    XIAD(4,2) = +1.d0
    XIAD(4,3) = 0.d0
!
    XIJM(1,1) = 0.d0
    XIJM(1,2) = -1.d0
    XIJM(1,3) = -1.d0
    XIJM(2,1) = 0.d0
    XIJM(2,2) = +1.d0
    XIJM(2,3) = -1.d0
    XIJM(3,1) = 0.d0
    XIJM(3,2) = +1.d0
    XIJM(3,3) = +1.d0
    XIJM(4,1) = 0.d0
    XIJM(4,2) = -1.d0
    XIJM(4,3) = +1.d0
!
! - Compute the decomposed derivatives of the jacobian in the covariant base
!
    call sshJacoDecoSB9(nb_node, geom ,&
                        J10    , J1ETA, J1ZETA, J1ETAZETA,&
                        J20    , J2XI , J2ZETA, J2XIZETA ,&
                        J30    , J3ETA, J3XI  , J3XIETA)
!
! - Bc0
!
    do I = 1, nb_node
        aux33 = 0.d0
        aux23 = 0.d0
        aux13 = 0.d0
        do AD = 1, 4
            call sshJacoColSB9(XIAD(AD,1:3),&
                               J10         , J1ETA, J1ZETA, J1ETAZETA,&
                               J20         , J2XI , J2ZETA, J2XIZETA ,&
                               J30         , J3ETA, J3XI  , J3XIETA  ,&
                               J1          , J2   , J3    )
            call sshShapeDeri(nb_node, XIAD(AD,1:3), dN_dXsi, dN_dEta, dN_dZeta)
            aux33 = aux33+dN_dZeta(I)*J3/4.d0
        enddo
        do EH = 1, 4
            call sshJacoColSB9(XIEH(EH,1:3),&
                               J10         , J1ETA, J1ZETA, J1ETAZETA,&
                               J20         , J2XI , J2ZETA, J2XIZETA ,&
                               J30         , J3ETA, J3XI  , J3XIETA  ,&
                               J1          , J2   , J3    )
            call sshShapeDeri(nb_node, XIEH(EH,1:3), dN_dXsi, dN_dEta, dN_dZeta)
            aux23 = aux23+(dN_dEta(I)*J3+dN_dZeta(I)*J2)/4.d0
        enddo
        do JM = 1, 4
            call sshJacoColSB9(XIJM(JM,1:3),&
                               J10         , J1ETA, J1ZETA, J1ETAZETA,&
                               J20         , J2XI , J2ZETA, J2XIZETA ,&
                               J30         , J3ETA, J3XI  , J3XIETA  ,&
                               J1          , J2   , J3    )
            call sshShapeDeri(nb_node, XIJM(JM,1:3), dN_dXsi, dN_dEta, dN_dZeta)
            aux13 = aux13+(dN_dXsi(I)*J3+dN_dZeta(I)*J1)/4.d0
        enddo
        C1=3*(I-1)+1
        C2=3*(I-1)+3
        Bc0(1,C1:C2) = g1(I)*J10
        Bc0(2,C1:C2) = g2(I)*J20
        Bc0(3,C1:C2) = aux33
        Bc0(4,C1:C2) = g1(I)*J20+g2(I)*J10
        Bc0(5,C1:C2) = aux13
        Bc0(6,C1:C2) = aux23
    enddo
!
! - BcZETA
!
    do I = 1, nb_node
        aux23 = 0.d0
        aux13 = 0.d0
        do EH = 1, 4
            call sshJacoColSB9(XIEH(EH,1:3),&
                               J10         , J1ETA, J1ZETA, J1ETAZETA,&
                               J20         , J2XI , J2ZETA, J2XIZETA ,&
                               J30         , J3ETA, J3XI  , J3XIETA  ,&
                               J1          , J2   , J3    )
            call sshShapeDeri(nb_node, XIEH(EH,1:3), dN_dXsi, dN_dEta, dN_dZeta)
            aux23 = aux23+XIEH(EH,3)*(dN_dEta(I)*J3+dN_dZeta(I)*J2)/4.d0
        enddo
        do JM = 1, 4
            call sshJacoColSB9(XIJM(JM,1:3),&
                               J10         , J1ETA, J1ZETA, J1ETAZETA,&
                               J20         , J2XI , J2ZETA, J2XIZETA ,&
                               J30         , J3ETA, J3XI  , J3XIETA  ,&
                               J1          , J2   , J3    )
            call sshShapeDeri(nb_node, XIJM(JM,1:3), dN_dXsi, dN_dEta, dN_dZeta)
            aux13 = aux13+XIJM(JM,3)*(dN_dXsi(I)*J3+dN_dZeta(I)*J1)/4.d0
        enddo
        C1=3*(I-1)+1
        C2=3*(I-1)+3
        aux12=h2(I)*J10+h3(I)*J20+g2(I)*J1ZETA+g1(I)*J2ZETA
        BcZETA(1,C1:C2) = h3(I)*J10+g1(I)*J1ZETA
        BcZETA(2,C1:C2) = h2(I)*J20+g2(I)*J2ZETA
        BcZETA(3,C1:C2) = 0.d0
        BcZETA(4,C1:C2) = aux12
        BcZETA(5,C1:C2) = aux13
        BcZETA(6,C1:C2) = aux23
    enddo
!
! - BcZETAZETA
!
    do I = 1, nb_node
        C1=3*(I-1)+1
        C2=3*(I-1)+3
        BcZETAZETA(1,C1:C2) = h3(I)*J1ZETA
        BcZETAZETA(2,C1:C2) = h2(I)*J2ZETA
        BcZETAZETA(3,C1:C2) = 0.d0
        BcZETAZETA(4,C1:C2) = h2(I)*J1ZETA+h3(I)*J2ZETA
        BcZETAZETA(5,C1:C2) = 0.d0
        BcZETAZETA(6,C1:C2) = 0.d0
    enddo
!
! - BcXI
!
    do I = 1, nb_node
        aux33 = 0.d0
        aux23 = 0.d0
        do AD = 1, 4
            call sshJacoColSB9(XIAD(AD,1:3),&
                               J10         , J1ETA, J1ZETA, J1ETAZETA,&
                               J20         , J2XI , J2ZETA, J2XIZETA ,&
                               J30         , J3ETA, J3XI  , J3XIETA  ,&
                               J1          , J2   , J3    )
            call sshShapeDeri(nb_node, XIAD(AD,1:3),dN_dXsi,dN_dEta,dN_dZeta)
            aux33 = aux33+XIAD(AD,1)*(dN_dZeta(I)*J3)/4.d0
        enddo
        do EH = 1, 4
            call sshJacoColSB9(XIEH(EH,1:3),&
                               J10         , J1ETA, J1ZETA, J1ETAZETA,&
                               J20         , J2XI , J2ZETA, J2XIZETA ,&
                               J30         , J3ETA, J3XI  , J3XIETA  ,&
                               J1          , J2   , J3    )
            call sshShapeDeri(nb_node, XIEH(EH,1:3),dN_dXsi,dN_dEta,dN_dZeta)
            aux23 = aux23+XIEH(EH,1)*(dN_dEta(I)*J3+dN_dZeta(I)*J2)/4.d0
        enddo
        C1=3*(I-1)+1
        C2=3*(I-1)+3
        BcXI(1,C1:C2) = 0.d0
        BcXI(2,C1:C2) = h1(I)*J20+g2(I)*J2XI
        BcXI(3,C1:C2) = aux33
        BcXI(4,C1:C2) = h1(I)*J10+g1(I)*J2XI
        BcXI(5,C1:C2) = 0.d0
        BcXI(6,C1:C2) = aux23
    enddo
!
! - BcETA
!
    do I = 1, nb_node
        aux33 = 0.d0
        aux13 = 0.d0
        do AD = 1, 4
            call sshJacoColSB9(XIAD(AD,1:3),&
                               J10         , J1ETA, J1ZETA, J1ETAZETA,&
                               J20         , J2XI , J2ZETA, J2XIZETA ,&
                               J30         , J3ETA, J3XI  , J3XIETA  ,&
                               J1          , J2   , J3    )
           call sshShapeDeri(nb_node, XIAD(AD,1:3),dN_dXsi,dN_dEta,dN_dZeta)
           aux33 = aux33+XIAD(AD,2)*(dN_dZeta(I)*J3)/4.d0
        enddo
        do JM = 1, 4
            call sshJacoColSB9(XIJM(JM,1:3),&
                               J10         , J1ETA, J1ZETA, J1ETAZETA,&
                               J20         , J2XI , J2ZETA, J2XIZETA ,&
                               J30         , J3ETA, J3XI  , J3XIETA  ,&
                               J1          , J2   , J3    )
            call sshShapeDeri(nb_node, XIJM(JM,1:3),dN_dXsi,dN_dEta,dN_dZeta)
            aux13 = aux13+XIJM(JM,2)*(dN_dXsi(I)*J3+dN_dZeta(I)*J1)/4.d0
        enddo
        C1=3*(I-1)+1
        C2=3*(I-1)+3
        BcETA(1,C1:C2) = h1(I)*J10+g1(I)*J1ETA
        BcETA(2,C1:C2) = 0.d0
        BcETA(3,C1:C2) = aux33
        BcETA(4,C1:C2) = h1(I)*J20+g2(I)*J1ETA
        BcETA(5,C1:C2) = aux13
        BcETA(6,C1:C2) = 0.d0
    enddo
!
! - BcETAZETA
!
    do I = 1, nb_node
        aux11 = h4(I)*J10+h3(I)*J1ETA+h1(I)*J1ZETA+g1(I)*J1ETAZETA
        aux12 = h4(I)*J20+h2(I)*J1ETA+h1(I)*J2ZETA+g2(I)*J1ETAZETA
        aux13 = 0.d0
        do JM = 1, 4
            call sshJacoColSB9(XIJM(JM,1:3),&
                               J10         , J1ETA, J1ZETA, J1ETAZETA,&
                               J20         , J2XI , J2ZETA, J2XIZETA ,&
                               J30         , J3ETA, J3XI  , J3XIETA  ,&
                               J1          , J2   , J3    )
            call sshShapeDeri(nb_node, XIJM(JM,1:3),dN_dXsi,dN_dEta,dN_dZeta)
            aux13 = aux13+XIJM(JM,2)*XIJM(JM,3)*(dN_dXsi(I)*J3+dN_dZeta(I)*J1)/4.d0
        enddo
        C1=3*(I-1)+1
        C2=3*(I-1)+3
        BcETAZETA(1,C1:C2) = aux11
        BcETAZETA(2,C1:C2) = 0.d0
        BcETAZETA(3,C1:C2) = 0.d0
        BcETAZETA(4,C1:C2) = aux12
        BcETAZETA(5,C1:C2) = aux13
        BcETAZETA(6,C1:C2) = 0.d0
    enddo
!
! - BcXIZETA
!
    do I = 1, nb_node
        aux22 = h4(I)*J20+h2(I)*J2XI+h1(I)*J2ZETA+g2(I)*J2XIZETA
        aux12 = h4(I)*J10+h3(I)*J2XI+h1(I)*J1ZETA+g1(I)*J2XIZETA
        aux23 = 0.d0
        do EH = 1, 4
            call sshJacoColSB9(XIEH(EH,1:3),&
                               J10         , J1ETA, J1ZETA, J1ETAZETA,&
                               J20         , J2XI , J2ZETA, J2XIZETA ,&
                               J30         , J3ETA, J3XI  , J3XIETA  ,&
                               J1          , J2   , J3    )
            call sshShapeDeri(nb_node, XIEH(EH,1:3),dN_dXsi,dN_dEta,dN_dZeta)
            aux23 = aux23+XIEH(EH,1)*XIEH(EH,3)*(dN_dEta(I)*J3+dN_dZeta(I)*J2)/4.d0
        enddo
        C1=3*(I-1)+1
        C2=3*(I-1)+3
        BcXIZETA(1,C1:C2) = 0.d0
        BcXIZETA(2,C1:C2) = aux22
        BcXIZETA(3,C1:C2) = 0.d0
        BcXIZETA(4,C1:C2) = aux12
        BcXIZETA(5,C1:C2) = 0.d0
        BcXIZETA(6,C1:C2) = aux23
    enddo
!
    !WRITE(6,*) 'Compute gradient matrix in covariant basis'
    !WRITE(6,*) ' * Bc0        :', bc0
    !WRITE(6,*) ' * BcXI       :', BcXI
    !WRITE(6,*) ' * BcETA      :', BcETA
    !WRITE(6,*) ' * BcZETA     :', BcZETA
    !WRITE(6,*) ' * BcXIZETA   :', BcXIZETA
    !WRITE(6,*) ' * BcETAZETA  :', BcETAZETA
    !WRITE(6,*) ' * BcZETAZETA :', BcZETAZETA
!
end subroutine
