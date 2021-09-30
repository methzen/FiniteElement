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
subroutine sshRigiGeomPtSB9(nbsig     , nb_node,&
                            nb_dof    , &
                            TZETA     , T      ,&
                            zeta      , sigm   ,&
                            matrGeomPt)
!
implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/sshShapeDeri.h"
#include "asterfort/sshTMatrSB9.h"
#include "asterfort/sshTMatrDecoSB9.h"
#include "asterfort/elrefe_info.h"
!
integer, intent(in) :: nb_node, nb_dof, nbsig
real(kind=8), intent(in) :: TZETA(6,6), T(6,6)
real(kind=8), intent(in) :: zeta, sigm(nbsig)
real(kind=8), intent(out) :: matrGeomPt(nb_dof, nb_dof)
!
! --------------------------------------------------------------------------------------------------
!
! Solid-shell element - SB9
!
! Compute geometric matrix at Gauss point
!
! --------------------------------------------------------------------------------------------------
!
! In  nbsig            : number of components o the stress tensor
! In  nb_node          : number of nodes of element without pinch node(s)
! In  nb_dof           : number of dof
! In  TZETA            : zeta components for T matrix
! In  T                : matrix relating the covariant and cartesian frames
! In  zeta             : zeta (out-of-plane parametric coordinate) of current Gauss point
! In  sigm             : stresses of current Gauss point
! Out matrGeomPt       : geometry rigidity matrix at current Gauss point
!
! --------------------------------------------------------------------------------------------------
!
    integer :: i, j, L1, L2, C1, C2, ad
    real(kind=8) :: XIEH(4,3), XIAD(4,3), XIJM(4,3)
    real(kind=8) :: dN_dXsi(8), dN_dEta(8), dN_dZeta(8)
    real(kind=8) :: G0(6), GZETA(6), GZETAZETA(6), const(6)
    real(kind=8) :: g1(8), g2(8), g3(8)
    real(kind=8) :: h1(8), h2(8), h3(8), h4(8)
    real(kind=8) :: aux33, aux23, aux13
    real(kind=8) :: auxz23, auxz13, auxz12
    real(kind=8) :: G9zeta(6), G9zzeta(6), G9zz(6)
    real(kind=8), parameter :: eye(3,3) = reshape((/1.d0, 0.d0, 0.d0,&
                                                    0.d0, 1.d0, 0.d0,&
                                                    0.d0, 0.d0, 1.d0/),(/3,3/))
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(nb_node .eq. 8)
    ASSERT(nb_dof .eq. 25)
    ASSERT(nbsig .eq. 6)
    matrGeomPt = 0.d0
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
    XIEH(2,2) =  0.d0
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
! - Compute geometric matrix
!
    do i = 1, nb_node
        do j = 1, nb_node
            L1 = 3*(i-1)+1
            L2 = 3*(i-1)+3
            C1 = 3*(j-1)+1
            C2 = 3*(j-1)+3
            aux33  = 0.d0
            aux23  = 0.d0
            aux13  = 0.d0
            auxz23 = 0.d0
            auxz13 = 0.d0
            do AD = 1, 4
                call sshShapeDeri(nb_node, XIAD(AD,1:3), dN_dXsi, dN_dEta, dN_dZeta)
                aux33  = aux33+dN_dZeta(i)*dN_dZeta(j)/4.d0
                call sshShapeDeri(nb_node, XIEH(AD,1:3), dN_dXsi, dN_dEta, dN_dZeta)
                aux23  = aux23+(dN_dEta(i)*dN_dZeta(j)+dN_dZeta(i)*dN_dEta(j))/4.d0
                auxz23 = auxz23+XIEH(AD,3)*(dN_dEta(i)*dN_dZeta(j)+dN_dZeta(i)*dN_dEta(j))/4.d0
                call sshShapeDeri(nb_node, XIJM(AD,1:3), dN_dXsi, dN_dEta, dN_dZeta)
                aux13  = aux13+(dN_dXsi(i)*dN_dZEta(j)+dN_dZeta(i)*dN_dXsi(j))/4.d0
                auxz13 = auxz13+XIJM(AD,3)*(dN_dXsi(i)*dN_dZeta(j)+dN_dZeta(i)*dN_dXsi(j))/4.d0
            end do
! --------- G0
            G0(1) = g1(i)*g1(j)
            G0(2) = g2(i)*g2(j)
            G0(3) = aux33
            G0(4) = g1(i)*g2(j)+g2(i)*g1(j)
            G0(5) = aux13
            G0(6) = aux23
! --------- GZETA
            auxz12   = h2(i)*g1(j)+h3(i)*g2(j)+g2(i)*h3(j)+g1(i)*h2(j)
            GZETA(1) = h3(i)*g1(j)+g1(i)*h3(j)
            GZETA(2) = h2(i)*g2(j)+g2(i)*h2(j)
            GZETA(3) = 0.d0
            GZETA(4) = auxz12
            GZETA(5) = auxz13
            GZETA(6) = auxz23
! --------- GZETAZETA
            GZETAZETA(1) = h3(i)*h3(j)
            GZETAZETA(2) = h2(i)*h2(j)
            GZETAZETA(3) = 0.d0
            GZETAZETA(4) = h2(i)*h3(j)+h2(j)*h3(i)
            GZETAZETA(5) = 0.d0
            GZETAZETA(6) = 0.d0
            const = matmul(T,G0)+&
                    zeta*(matmul(T,GZETA) + matmul(TZETA,G0))+&
                    zeta*zeta*(matmul(T,GZETAZETA) + matmul(TZETA,GZETA))
            matrGeomPt(L1:L2,C1:C2) = sum(const*sigm)*eye
        enddo
    enddo
!
! - Pinch node
!
    G9zeta  = 0.d0
    G9zzeta = 0.d0
    G9zz    = 0.d0
    do i = 1, 8
        G9zeta(3)  = -2.d0*g3(i)
        G9zeta(5)  = -2.d0*g1(i)
        G9zeta(6)  = -2.d0*g2(i)
        G9zzeta(5) = -2.d0*h3(i)
        G9zzeta(6) = -2.d0*h2(i)
        const = zeta*matmul(T,G9zeta)+&
                zeta*zeta*(matmul(T,G9zzeta) + matmul(TZETA,G9zeta))
        matrGeomPt(25,3*(i-1)+3) = sum(const*sigm)
        matrGeomPt(3*(i-1)+3,25) = sum(const*sigm)
    enddo
    G9zz(3) = 4.d0
    const   = matmul(T,G9zz)*zeta*zeta
    matrGeomPt(3*nb_node+1, 3*nb_node+1) = sum(const*sigm)
!
end subroutine
