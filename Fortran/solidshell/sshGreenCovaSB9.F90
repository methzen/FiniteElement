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
subroutine sshGreenCovaSB9(nb_node  , nb_dof ,&
                           geomInit , disp   ,&
                           Ec0_     , EcZETA_, EcZETAZETA_,&
                           EcXI_    , EcETA_ , EcETAZETA_ ,&
                           EcXIZETA_)
!
implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/sshJacoDecoSB9.h"
#include "asterfort/sshDispDeriCovaSB9.h"
#include "asterfort/sshJacoColSB9.h"
#include "asterfort/sshDispDeriCovaColSB9.h"
!
integer, intent(in) :: nb_node, nb_dof
real(kind=8), intent(in) :: geomInit(3*nb_node), disp(nb_dof)
real(kind=8), optional, intent(out) :: Ec0_(6), EcZETA_(6), EcZETAZETA_(6)
real(kind=8), optional, intent(out) :: EcXI_(6), EcETA_(6), EcETAZETA_(6), EcXIZETA_(6)
!
! --------------------------------------------------------------------------------------------------
!
! Solid-shell element - SB9
!
! Compute Green-Lagrange strains in covariant basis
!
! --------------------------------------------------------------------------------------------------
!
! Green-Lagrange strains in covariant frame
!
!   Ec  = Ec0 + xi.EcXI + eta.EcETA + zeta.EcZETA +
!         xi.eta.EcXIETA + eta.zeta.EcETAZETA + xi.zeta.EcXIZETA
!
! In  nb_node          : number of nodes of element without pinch node(s)
! In  nb_dof           : number of dof
! In  geomInit         : initial coordinates of element
! In  disp             : displacements of element
! Out Ec0              : constant part of Green-Lagrange strains in covariant frame
! Out EcXI             : xi part of Green-Lagrange strains in covariant frame
! Out EcETA            : eta part of Green-Lagrange strains in covariant frame
! Out EcZETA           : zeta part of Green-Lagrange strains in covariant frame
! Out EcXIETA          : (xi,eta) part of Green-Lagrange strains in covariant frame
! Out EcETAZETA        : (eta,zeta) part of Green-Lagrange strains in covariant frame
! Out EcXIZETA         : (xi,zeta) part of Green-Lagrange strains in covariant frame
!
! --------------------------------------------------------------------------------------------------
!
    integer :: AD, EH, JM
    real(kind=8) :: XIEH(4,3), XIAD(4,3) ,XIJM(4,3)
    real(kind=8) :: J10(3), J1ETA(3), J1ZETA(3), J1ETAZETA(3)
    real(kind=8) :: J20(3), J2XI(3), J2ZETA(3), J2XIZETA(3)
    real(kind=8) :: J30(3), J3ETA(3), J3XI(3), J3XIETA(3)
    real(kind=8) :: D10(3), D1ETA(3), D1ZETA(3), D1ETAZETA(3)
    real(kind=8) :: D20(3), D2XI(3), D2ZETA(3), D2XIZETA(3)
    real(kind=8) :: D30(3), D3ETA(3), D3XI(3), D3XIETA(3)
    real(kind=8) :: J1(3), J2(3), J3(3), D1(3), D2(3), D3(3)
    real(kind=8) :: EcXI(6), EcETA(6), EcETAZETA(6), EcXIZETA(6)
    real(kind=8) :: Ec0(6), EcZETA(6), EcZETAZETA(6)
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(nb_node .eq. 8)
    ASSERT(nb_dof .eq. 25)
!
    Ec0        = 0.d0
    EcZETA     = 0.d0
    EcZETAZETA = 0.d0
    EcXI       = 0.d0
    EcETA      = 0.d0
    EcETAZETA  = 0.d0
    EcXIZETA   = 0.d0
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
    call sshJacoDecoSB9(nb_node, geomInit,&
                                 J10    , J1ETA    , J1ZETA, J1ETAZETA,&
                                 J20    , J2XI     , J2ZETA, J2XIZETA ,&
                                 J30    , J3ETA    , J3XI  , J3XIETA)
!
! - Compute the decomposed derivatives of the displacements in the covariant base
!
    call sshDispDeriCovaSB9(nb_node, nb_dof, disp  ,&
                            D10    , D1ETA , D1ZETA, D1ETAZETA,&
                            D20    , D2XI  , D2ZETA, D2XIZETA ,&
                            D30    , D3ETA , D3XI  , D3XIETA  )
!
! - Compute Ec0
!
    Ec0(1) = sum(J10*D10)+0.5d0*sum(D10*D10)
    Ec0(2) = sum(J20*D20)+0.5d0*sum(D20*D20)
    Ec0(4) = sum(J10*D20)+sum(J20*D10)+sum(D10*D20)
!
! - Compute EcZETA
!
    EcZETA(1) = sum(J10*D1ZETA)+sum(D10*J1ZETA)+sum(D10*D1ZETA)
    EcZETA(2) = sum(J20*D2ZETA)+sum(D20*J2ZETA)+sum(D20*D2ZETA)
    EcZETA(4) = sum(J10*D2ZETA)+sum(D20*J1ZETA)+sum(J20*D1ZETA)+&
                sum(D10*J2ZETA)+sum(D10*D2ZETA)+sum(D20*D1ZETA)
!
! - Compute EcZETAZETA
!
    EcZETAZETA(1) = sum(J1ZETA*D1ZETA)+0.5d0*sum(D1ZETA*D1ZETA)
    EcZETAZETA(2) = sum(J2ZETA*D2ZETA)+0.5d0*sum(D2ZETA*D2ZETA)
    EcZETAZETA(4) = sum(J1ZETA*D2ZETA)+sum(J2ZETA*D1ZETA)+sum(D2ZETA*D1ZETA)
!
! - Compute EcXI (stabilization)
!
    EcXI(2) = sum(J20*D2XI)+sum(J2XI*D20)+sum(D2XI*D20)
    EcXI(4) = sum(J10*D2XI)+sum(J2XI*D10)+sum(D2XI*D10)
!
! - Compute EcETA (stabilization)
!
    EcETA(1) = sum(J10*D1ETA)+sum(D10*J1ETA)+sum(D10*D1ETA)
    EcETA(4) = sum(J20*D1ETA)+sum(D20*J1ETA)+sum(D20*D1ETA)
!
! - Compute ETAZETA (stabilization)
!
    EcETAZETA(1) = sum(J10*D1ETAZETA)+sum(J1ETA*D1ZETA)+sum(J1ZETA*D1ETA)+&
                   sum(D10*J1ETAZETA)+sum(D10*D1ETAZETA)+sum(D1ETA*D1ZETA)
    EcETAZETA(4) = sum(J20*D1ETAZETA)+sum(J1ETA*D2ZETA)+sum(J2ZETA*D1ETA)+&
                   sum(D20*J1ETAZETA)+sum(D20*D1ETAZETA)+sum(D1ETA*D2ZETA)
!
! - Compute XIZETA (stabilization)
!
    EcXIZETA(2) = sum(J20*D2XIZETA)+sum(J2XI*D2ZETA)+sum(J2ZETA*D2XI)+&
                  sum(D20*J2XIZETA)+sum(D20*D2XIZETA)+sum(D2XI*D2ZETA)
    EcXIZETA(4) = sum(J10*D2XIZETA)+sum(J2XI*D1ZETA)+sum(J1ZETA*D2XI)+&
                  sum(D10*J2XIZETA)+sum(D10*D2XIZETA)+sum(D2XI*D1ZETA)
!
! - Assumed strains
!
    do AD = 1, 4
        call sshJacoColSB9(XIAD(AD,1:3),&
                                   J10, J1ETA, J1ZETA, J1ETAZETA,&
                                   J20, J2XI , J2ZETA, J2XIZETA ,&
                                   J30, J3ETA, J3XI  , J3XIETA  ,&
                                   J1 , J2   , J3    )
        call sshDispDeriCovaColSB9(XIAD(AD,1:3),&
                                   D10         , D1ETA, D1ZETA, D1ETAZETA,&
                                   D20         , D2XI , D2ZETA, D2XIZETA ,&
                                   D30         , D3ETA, D3XI  , D3XIETA  ,&
                                   D1          , D2   , D3    )
        Ec0(3)   = Ec0(3)+(sum(J3*D3)+0.5d0*sum(D3*D3))/4.d0
        EcXI(3)  = EcXI(3)+XIAD(AD,1)*(sum(J3*D3)+0.5d0*sum(D3*D3))/4.d0
        EcETA(3) = EcETA(3)+XIAD(AD,2)*(sum(J3*D3)+0.5d0*sum(D3*D3))/4.d0
    enddo
!
    do EH = 1, 4
        call sshJacoColSB9(XIEH(EH,1:3),&
                                   J10, J1ETA, J1ZETA, J1ETAZETA,&
                                   J20, J2XI , J2ZETA, J2XIZETA ,&
                                   J30, J3ETA, J3XI  , J3XIETA  ,&
                                   J1 , J2   , J3    )
        call sshDispDeriCovaColSB9(XIEH(EH,1:3),&
                                   D10         , D1ETA, D1ZETA, D1ETAZETA,&
                                   D20         , D2XI , D2ZETA, D2XIZETA ,&
                                   D30         , D3ETA, D3XI  , D3XIETA  ,&
                                   D1          , D2   , D3    )
        Ec0(6)      = Ec0(6)+(sum(J3*D2)+sum(J2*D3)+sum(D3*D2))/4.d0
        EcZETA(6)   = EcZETA(6)+XIEH(EH,3)*(sum(J3*D2)+sum(J2*D3)+sum(D3*D2))/4.d0
        EcXI(6)     = EcXI(6)+XIEH(EH,1)*(sum(J2*D3)+sum(J3*D2)+sum(D2*D3))/4.d0
        EcXIZETA(6) = EcXIZETA(6)+XIEH(EH,1)*XIEH(EH,3)*(sum(J2*D3)+sum(J3*D2)+sum(D2*D3))/4.d0
    enddo
!
    do JM = 1, 4
        call sshJacoColSB9(XIJM(JM,1:3),&
                                   J10, J1ETA, J1ZETA, J1ETAZETA,&
                                   J20, J2XI , J2ZETA, J2XIZETA ,&
                                   J30, J3ETA, J3XI  , J3XIETA  ,&
                                   J1 , J2   , J3    )
        call sshDispDeriCovaColSB9(XIJM(JM,1:3),&
                                   D10         , D1ETA, D1ZETA, D1ETAZETA,&
                                   D20         , D2XI , D2ZETA, D2XIZETA ,&
                                   D30         , D3ETA, D3XI  , D3XIETA  ,&
                                   D1          , D2   , D3    )
       Ec0(5)       = Ec0(5)+(sum(J3*D1)+sum(J1*D3)+sum(D3*D1))/4.d0
       EcZETA(5)    = EcZETA(5)+XIJM(JM,3)*(sum(J3*D1)+sum(J1*D3)+sum(D3*D1))/4.d0
       EcETA(5)     = EcETA(5)+XIJM(JM,2)*(sum(J1*D3)+sum(J3*D1)+sum(D1*D3))/4.d0
       EcETAZETA(5) = EcETAZETA(5)+XIJM(JM,2)*XIJM(JM,3)*(sum(J1*D3)+sum(J3*D1)+sum(D1*D3))/4.d0
    enddo
!
    if (present(EcXI_)) then
        EcXI_ = EcXI
    endif
    if (present(EcETA_)) then
        EcETA_ = EcETA
    endif
    if (present(EcETAZETA_)) then
        EcETAZETA_ = EcETAZETA
    endif
    if (present(EcXIZETA_)) then
        EcXIZETA_ = EcXIZETA
    endif
    if (present(Ec0_)) then
        Ec0_ = Ec0
    endif
    if (present(EcZETA_)) then
        EcZETA_ = EcZETA
    endif
    if (present(EcZETAZETA_)) then
        EcZETAZETA_ = EcZETAZETA
    endif

!
    !WRITE(6,*) 'Déformation de Green-Lagrange (covariantes): ',Ec0      , EcZETA, EcZETAZETA,&
    !                                                           EcXI     , EcETA , EcETAZETA ,&
    !                                                           EcXIZETA
!
end subroutine
