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
subroutine sshNLStabMatrGeomSB9(nb_node    , nb_dof,&
                                geomInit   , disp      , Ueff,&
                                BXIdev     , BETAdev   ,&
                                BETAZETAdev, BXIZETAdev,&
                                Kgstab)
!
implicit none
!
#include "asterf_types.h"
#include "asterfort/sshTMatrDecoSB9.h"
#include "asterfort/sshTMatrSB9.h"
#include "asterfort/sshShapeDeri.h"
#include "asterfort/assert.h"
#include "asterfort/sshNLStabSigmSB9.h"
!
integer, intent(in) :: nb_node, nb_dof
real(kind=8), intent(in) :: geomInit(3*nb_node), disp(nb_dof), Ueff
real(kind=8), intent(in) :: BXIdev(6,3*nb_node), BETAdev(6,3*nb_node)
real(kind=8), intent(in) :: BXIZETAdev(6,3*nb_node), BETAZETAdev(6,3*nb_node)
real(kind=8), intent(out) :: Kgstab(3*nb_node,3*nb_node)
!
! --------------------------------------------------------------------------------------------------
!
! Solid-shell element - SB9
!
! Compute stabilization for geometric rigidity
!
! --------------------------------------------------------------------------------------------------
!
! In  nb_node          : number of nodes of element without pinch node(s)
! In  nb_dof           : number of dof
! In  geomInit        : initial coordinates of element
! In  disp             : displacements of element
! In  Ueff             : effective shear modulus
! In  BXIdev           : xi part of gradient matrix (deviatoric part)
! In  BETAdev          : eta part of gradient matrix (deviatoric part)
! In  BETAZETAdev      : eta x zeta part of gradient matrix (deviatoric part)
! In  BXIZETAdev       : xi x zeta part of gradient matrix (deviatoric part)
! Out Kgstab           : stabilisation for geometric rigidity matrix
!
! --------------------------------------------------------------------------------------------------
!
    integer :: i, j, L1, L2, C1, C2, AD, EH, JM
    real(kind=8) :: XIEH(4,3), XIAD(4,3), XIJM(4,3)
    real(kind=8) :: dN_dXsi(8), dN_dEta(8), dN_dZeta(8)
    real(kind=8) :: aux13, aux23, aux33, aux12, aux11, aux22
    real(kind=8) :: SHXI(6), SHETA(6), SHETAZETA(6), SHXIZETA(6)
    real(kind=8) :: TXI(6,6), TETA(6,6), TZETA(6,6), T(6,6)
    real(kind=8) :: const
    real(kind=8)::  g1(8), g2(8), g3(8)
    real(kind=8)::  h1(8), h2(8), h3(8), h4(8)
    real(kind=8) :: Gc0(6), Gczeta(6), GcXI(6)
    real(kind=8) :: GcETA(6), GcETAZETA(6), GcXIZETA(6)
    real(kind=8) :: GXI(6), GETA(6), GETAZETA(6), GXIZETA(6)
    aster_logical, parameter :: lLarge = ASTER_TRUE
    real(kind=8), parameter :: eye(3,3) = reshape((/1.d0, 0.d0, 0.d0,&
                                                    0.d0, 1.d0, 0.d0,&
                                                    0.d0, 0.d0, 1.d0/),(/3,3/))
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(nb_node .eq. 8)
    ASSERT(nb_dof .eq. 25)
    Kgstab = 0.d0
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
    XIEH(2,2) = +0.d0
    XIEH(2,3) = -1.d0
    XIEH(3,1) = +1.d0
    XIEH(3,2) =  0.d0
    XIEH(3,3) = +1.d0
    XIEH(4,1) = -1.d0
    XIEH(4,2) =  0.d0
    XIEH(4,3) = +1.d0
!
    XIAD(1,1) = -1.d0
    XIAD(1,2) = -1.d0
    XIAD(1,3) =  0.d0
    XIAD(2,1) = +1.d0
    XIAD(2,2) = -1.d0
    XIAD(2,3) =  0.d0
    XIAD(3,1) = +1.d0
    XIAD(3,2) = +1.d0
    XIAD(3,3) =  0.d0
    XIAD(4,1) = -1.d0
    XIAD(4,2) = +1.d0
    XIAD(4,3) =  0.d0
!
    XIJM(1,1) =  0.d0
    XIJM(1,2) = -1.d0
    XIJM(1,3) = -1.d0
    XIJM(2,1) =  0.d0
    XIJM(2,2) = +1.d0
    XIJM(2,3) = -1.d0
    XIJM(3,1) =  0.d0
    XIJM(3,2) = +1.d0
    XIJM(3,3) = +1.d0
    XIJM(4,1) =  0.d0
    XIJM(4,2) = -1.d0
    XIJM(4,3) = +1.d0
!
! - Compute the decomposition of the inverse Jacobian matrix
!
    call sshTMatrDecoSB9(nb_node, geomInit, TXI, TETA, TZETA)
!
! - Compute T matrix (matrix relating the covariant and cartesian frames) at center of element
!
    call sshTMatrSB9(nb_node, geomInit, (/0.d0,0.d0,0.d0/), T)
!
! - Compute stabilisation stresses
!
    call sshNLStabSigmSB9(lLarge     ,&
                          nb_node    , nb_dof    ,&
                          geomInit   , disp      , Ueff,&
                          BXIdev     , BETAdev   ,&
                          BETAZETAdev, BXIZETAdev,&
                          SHXI       , SHETA     ,&
                          SHETAZETA  , SHXIZETA)
!
! - Compute matrix
!
    do i = 1, nb_node
        do j = 1, nb_node
            L1 = 3*(i-1)+1
            L2 = 3*(i-1)+3
            C1 = 3*(j-1)+1
            C2 = 3*(j-1)+3
! --------- G0
            aux33 = 0.d0
            aux23 = 0.d0
            aux13 = 0.d0
            do AD = 1, 4
                call sshShapeDeri(nb_node, XIAD(AD,1:3), dN_dXsi, dN_dEta, dN_dZeta)
                aux33 = aux33+dN_dZeta(i)*dN_dZeta(j)
            end do
            do EH = 1, 4
                call sshShapeDeri(nb_node, XIEH(EH,1:3), dN_dXsi, dN_dEta, dN_dZeta)
                aux23 = aux23+(dN_dEta(i)*dN_dZeta(j)+dN_dZeta(i)*dN_dEta(j))/4.d0
            end do
            do JM = 1, 4
                call sshShapeDeri(nb_node, XIJM(JM,1:3), dN_dXsi, dN_dEta, dN_dZeta)
                aux13 = aux13+(dN_dXsi(i)*dN_dZeta(j)+dN_dZeta(i)*dN_dXsi(j))/4.d0
            end do
            Gc0(1) = g1(i)*g1(j)
            Gc0(2) = g2(i)*g2(j)
            Gc0(3) = aux33
            Gc0(4) = g1(i)*g2(j)+g2(i)*g1(j)
            Gc0(5) = aux13
            Gc0(6) = aux23
! --------- GZETA
            aux23 = 0.d0
            aux13 = 0.d0
            do EH = 1, 4
                call sshShapeDeri(nb_node, XIEH(EH,1:3), dN_dXsi, dN_dEta, dN_dZeta)
                aux23 = aux23+XIEH(EH,3)*(dN_dEta(i)*dN_dZeta(j)+dN_dZeta(i)*dN_dEta(j))/4.d0
            end do
            do JM = 1, 4
                call sshShapeDeri(nb_node, XIJM(JM,1:3), dN_dXsi, dN_dEta, dN_dZeta)
                aux13 = aux13+XIJM(JM,3)*(dN_dXsi(i)*dN_dZeta(j)+dN_dZeta(i)*dN_dXsi(j))/4.d0
            end do
            aux12     = h2(i)*g1(j)+h3(i)*g2(j)+g2(i)*h3(j)+g1(i)*h2(j)
            GcZETA(1) = h3(i)*g1(j)+g1(i)*h3(j)
            GcZETA(2) = h2(i)*g2(j)+g2(i)*h2(j)
            GcZETA(3) = 0.d0
            GcZETA(4) = aux12
            GcZETA(5) = aux13
            GcZETA(6) = aux23
! --------- GcXI
            aux33 = 0.d0
            aux23  =0.d0
            do AD = 1, 4
                call sshShapeDeri(nb_node, XIAD(AD,1:3), dN_dXsi, dN_dEta, dN_dZeta)
                aux33 = aux33+XIAD(AD,1)*(dN_dZeta(i)*dN_dZeta(j))/4.d0
            end do
            do EH = 1, 4
                call sshShapeDeri(nb_node, XIEH(EH,1:3),dN_dXsi,dN_dEta,dN_dZeta)
                aux23 = aux23+XIEH(EH,1)*(dN_dEta(i)*dN_dZeta(j)+dN_dZeta(i)*dN_dEta(j))/4.d0
            end do
            GcXI(1) = 0.d0
            GcXI(2) = h1(i)*g2(j)+g2(i)*h1(j)
            GcXI(3) = aux33
            GcXI(4) = h1(i)*g1(j)+g1(i)*h1(j)
            GcXI(5) = 0.d0
            GcXI(6) = aux23
! --------- Term GcETA
            aux33 = 0.d0
            aux13 = 0.d0
            do AD = 1, 4
                call sshShapeDeri(nb_node, XIAD(AD,1:3), dN_dXsi, dN_dEta, dN_dZeta)
                aux33 = aux33+XIAD(AD,2)*(dN_dZeta(i)*dN_dZeta(j))/4.d0
            end do
            do JM = 1, 4
                call sshShapeDeri(nb_node, XIJM(JM,1:3), dN_dXsi, dN_dEta, dN_dZeta)
                aux13 = aux13+XIJM(JM,2)*(dN_dXsi(i)*dN_dZeta(j)+dN_dZeta(i)*dN_dXsi(j))/4.d0
            end do
            GcETA(1) = h1(i)*g1(j)+g1(i)*h1(j)
            GcETA(2) = 0.d0
            GcETA(3) = aux33
            GcETA(4) = h1(i)*g2(j)+g2(i)*h1(j)
            GcETA(5) = aux13
            GcETA(6) = 0.d0
! --------- GcETAZETA
            aux11 = h4(i)*g1(j)+h3(i)*h1(j)+h1(i)*h3(j)+g1(i)*h4(j)
            aux12 = h4(i)*g2(j)+h2(i)*h1(j)+h1(i)*h2(j)+g2(i)*h4(j)
            aux13 = 0.d0
            do JM = 1, 4
                call sshShapeDeri(nb_node, XIJM(JM,1:3), dN_dXsi, dN_dEta, dN_dZeta)
                aux13 = aux13+&
                        XIJM(JM,2)*XIJM(JM,3)*(dN_dXsi(i)*dN_dZeta(j)+dN_dXsi(j)*dN_dZeta(i))/4.d0
            end do
            GcETAZETA(1) = aux11
            GcETAZETA(2) = 0.d0
            GcETAZETA(3) = 0.d0
            GcETAZETA(4) = aux12
            GcETAZETA(5) = aux13
            GcETAZETA(6) = 0.d0
! --------- GcXIZETA
            aux22 = h4(i)*g2(j)+h2(i)*h1(j)+h1(i)*h2(j)+g2(i)*h4(j)
            aux12 = h4(i)*g1(j)+h3(i)*h1(j)+h1(i)*h3(j)+g1(i)*h4(j)
            aux23 = 0.d0
            do EH = 1, 4
                call sshShapeDeri(nb_node, XIEH(EH,1:3), dN_dXsi, dN_dEta, dN_dZeta)
                aux23 = aux23+&
                        XIEH(EH,1)*XIEH(EH,3)*(dN_dEta(i)*dN_dZeta(j)+dN_dEta(j)*dN_dZeta(i))/4.d0
            end do
            GcXIZETA(1) = 0.d0
            GcXIZETA(2) = aux22
            GcXIZETA(3) = 0.d0
            GcXIZETA(4) = aux12
            GcXIZETA(5) = 0.d0
            GcXIZETA(6) = aux23
!
            GXI      = matmul(T,GcXI)+matmul(TXI,Gc0)
            GETA     = matmul(T,GcETA)+matmul(TETA,Gc0)
            GETAZETA = matmul(T,GcETAZETA)+matmul(TETA,GcZETA)+matmul(TZETA,GcETA)
            GXIZETA  = matmul(T,GcXIZETA)+matmul(TXI,GcZETA)+matmul(TZETA,GcXI)
            const = sum(SHXI*GXI)*(8.d0/3.d0)+sum(SHETA*GETA)*(8.d0/3.d0)+&
                   (sum(SHETAZETA*GETAZETA)+sum(SHXIZETA*GXIZETA))*(8.d0/9.d0)
            Kgstab(L1:L2,C1:C2) = const*eye
        end do
    end do
!
    !WRITE(6,*) 'Compute stabilisation for geometric rigidity matrix: ',Kgstab
!
end subroutine
