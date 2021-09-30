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
subroutine sshRigiMatrSB7(npg        , nb_node , nb_dof  ,&
                          jv_geom    , jv_coopg, jv_poids,&
                          hookeMatrix, matrRigi)
!
implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/sshBTSB.h"
#include "asterfort/sshGradMatrSB7.h"
!
integer, intent(in) :: npg, nb_node, nb_dof
integer, intent(in) :: jv_geom, jv_poids, jv_coopg
real(kind=8), intent(in) :: hookeMatrix(6, 6)
real(kind=8), intent(out) :: matrRigi(25, 25)
!
! --------------------------------------------------------------------------------------------------
!
! Solid-shell element - SB7
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
    integer :: i, j, kpg
    real(kind=8) :: geomInit(18)
    real(kind=8) :: zeta, poids, jacobian
    real(kind=8) :: B(6, 19), BStab(2, 18)
    real(kind=8) :: tBDB(19, 19), tBsDBs(18, 18)
    real(kind=8) :: Kstab(19, 19), Dsz(2, 2)
    real(kind=8), parameter :: eye(2,2) = reshape((/1.d0, 0.d0,&
                                                    0.d0, 1.d0/),(/2,2/))
    real(kind=8), parameter :: stabCoef = 1200.d0
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(nb_node .eq. 6)
    ASSERT(nb_dof .eq. 19)
!
    matrRigi = 0.d0
!
! - Get initial geometry
!
    do i = 1, 3*nb_node
        geomInit(i) = zr(jv_geom+i-1)
    end do
!
! - Loop on Gauss points
!
    do kpg = 1, npg
        zeta  = zr(jv_coopg-1+3*kpg)
        poids = zr(jv_poids+kpg-1)
! ----- Compute B matrix at current Gauss point
        call sshGradMatrSB7(nb_node, nb_dof, geomInit, zeta,&
                            B      , BStab , jacobian)
! ----- Compute rigidity term
        !WRITE(6,*) 'hookeMatrix: ',sum(hookeMatrix)
        call sshBTSB(hookeMatrix, 6, nb_dof, B, tBDB)
        !WRITE(6,*) 'tBDB: ',sum(tBDB)
! ----- Compute stabilization term
        Kstab = 0.d0
        Dsz   = 0.d0
        Dsz   = eye*hookeMatrix(3, 3)/stabCoef
        !WRITE(6,*) 'Dsz: ',sum(Dsz)
        call sshBTSB(Dsz, 2, 18, BStab, tBsDBs)
        !WRITE(6,*) 'tBsDBs: ',sum(tBsDBs)
        Kstab(1:18, 1:18) = tBsDBs
        !WRITE(6,*) 'Kstab: ',sum(Kstab)
! ----- Compute matrix
        do i = 1, nb_dof
            do j = 1, nb_dof
                matrRigi(i,j) = matrRigi(i,j) + jacobian*poids*(tBDB(i,j)+Kstab(i,j))
            end do
        end do
    end do
    !WRITE(6,*) 'matrRigi: ',sum(matrRigi)
!
end subroutine
