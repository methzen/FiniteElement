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
subroutine sshRigiMatrSB9(npg        , nb_node , nb_dof  ,&
                          jv_geom    , jv_coopg, jv_poids,&
                          hookeMatrix, matrRigi)
!
implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dfdm3d.h"
#include "asterfort/sshBTSB.h"
#include "asterfort/sshGradMatrCartSB9.h"
#include "asterfort/sshGradMatrSB9.h"
#include "asterfort/sshStabMatrSB9.h"
!
integer, intent(in) :: npg, nb_node, nb_dof
integer, intent(in) :: jv_geom, jv_poids, jv_coopg
real(kind=8), intent(in) :: hookeMatrix(6, 6)
real(kind=8), intent(out) :: matrRigi(nb_dof, nb_dof)
!
! --------------------------------------------------------------------------------------------------
!
! Solid-shell element - SB9
!
! Compute rigidity matrix
!
! --------------------------------------------------------------------------------------------------
!
! In  npg              : number of Gauss points
! In  nb_node          : number of nodes of element without pinch node(s)
! In  nb_dof           : number of dof
! In  jv_geom          : JEVEUX adress to geometry
! In  jv_poids         : JEVEUX adress to weight of Gauss points
! In  jv_coopg         : JEVEUX adress to coordinates of Gauss points
! In  hookeMatrix      : elastic matrix
! Out matrRigi         : rigidity matrix
!
! --------------------------------------------------------------------------------------------------
!
    integer :: kpg, i, j
    real(kind=8) :: zeta, poids, aux1, aux2, det
    real(kind=8) :: geomInit(24), disp(25), Ueff
    real(kind=8) :: tBDB(25,25)
    real(kind=8) :: B0(6,24), B(6,25)
    real(kind=8) :: SXI(24,24), SETA(24,24), SETAZETA(24,24), SXIZETA(24,24)
    real(kind=8) :: BXI(6,24), BETA(6,24), BZETA(6,24)
    real(kind=8) :: BXIZETA(6,24), BETAZETA(6,24), BZETAZETA(6,24)
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(nb_node .eq. 8)
    ASSERT(nb_dof .eq. 25)
!
    matrRigi = 0.d0
    disp     = 0.d0
!
! - Get initial geometry
!
    do i = 1, 3*nb_node
        geomInit(i) = zr(jv_geom+i-1)
    end do
!
! - Compute the deformation gradient in cartesian base
!
    call sshGradMatrCartSB9(nb_node  , geomInit,&
                            nb_dof   , disp    ,&
                            B0       , BZETA   , BZETAZETA,&
                            BXI      , BETA    , &
                            BXIZETA  , BETAZETA, &
                            det)
!
! - Loop on Gauss points
!
    do kpg = 1, npg
        zeta  = zr(jv_coopg-1+3*kpg)
        poids = zr(jv_poids+kpg-1)
! ----- Compute B matrix at current Gauss point
        call sshGradMatrSB9(nb_node, nb_dof   ,&
                            zeta   , geomInit,&
                            B0     , BZETA    , BZETAZETA,&
                            B)
! ----- Compute product tBSB
        call sshBTSB(hookeMatrix, 6, nb_dof, B, tBDB)
        do i = 1, nb_dof
            do j = 1, nb_dof
                matrRigi(i,j) = matrRigi(i,j) + det*poids*tBDB(i,j)
            enddo
        enddo
    end do
!
! - Compute stabilization matrix (elastic)
!
    Ueff = hookeMatrix(5,5)
    !print*,'Ueff',Ueff
    call sshStabMatrSB9(nb_node, Ueff    ,&
                        BXI    , BETA    ,&
                        BXIZETA, BETAZETA,&
                        SXI    , SETA    ,&
                        SXIZETA, SETAZETA)
!
! - Compute matrix
!
    aux1  = det*8.d0/3.d0
    aux2  = det*8.d0/9.d0
    do i = 1, 3*nb_node
        do j = 1, 3*nb_node
            matrRigi(i,j) = matrRigi(i,j)+&
                            (SXI(i,j)+SETA(i,j))*aux1+&
                            (SXIZETA(i,j)+SETAZETA(i,j))*aux2
        enddo
    enddo
!
end subroutine
