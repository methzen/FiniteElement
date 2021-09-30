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
subroutine sshSiefElgaSB7(npg     , nbsig   ,&
                          nb_node , nb_dof  ,&
                          j_mater , jv_coopg,&
                          geom    , disp    ,&
                          sigm)
!
implicit none
!
#include "jeveux.h"
#include "asterc/r8vide.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/dmat3d.h"
#include "asterfort/sshGradMatrSB7.h"
!
integer, intent(in) :: npg, nb_node, nb_dof, nbsig
integer, intent(in) :: jv_coopg, j_mater
real(kind=8), intent(in) :: geom(3*nb_node), disp(nb_dof)
real(kind=8), intent(out) :: sigm(nbsig*npg)
!
! --------------------------------------------------------------------------------------------------
!
! Solid-shell element - SB7
!
! Compute SIEF_ELGA
!
! --------------------------------------------------------------------------------------------------
!
! In  npg              : number of Gauss points
! In  nbsig            : number of components o the stress tensor
! In  nb_node          : number of nodes of element without pinch node(s)
! In  nb_dof           : number of dof
! In  j_mater          : JEVEUX address for material properties
! In  jv_coopg         : JEVEUX adress to coordinates of Gauss points
! In  geom             : coordinates of element
! In  disp             : displacements of element
! Out sigm             : stresses
!
! --------------------------------------------------------------------------------------------------
!
    integer :: kpg
    real(kind=8) :: hookeMatrix(6, 6)
    real(kind=8) :: repere(7), xyzgau(3), epsi(6)
    real(kind=8) :: B(6, 19), zeta
    character(len=4), parameter :: fami = 'RIGI'
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(nb_node .eq. 6)
    ASSERT(nb_dof .eq. 19)
    ASSERT(nbsig .eq. 6)
!
    hookeMatrix = 0.d0
    repere      = 0.d0
    xyzgau      = 0.d0
    sigm        = 0.d0
!
! - Get elastic matrix
!
    call dmat3d(fami, j_mater , r8vide(), '+', 1,&
                1   , repere, xyzgau, hookeMatrix)
!
! - Loop on Gauss points
!
    do kpg = 1, npg
        zeta  = zr(jv_coopg-1+3*kpg)
! ----- Compute B matrix
        call sshGradMatrSB7(nb_node, nb_dof, geom, zeta, B)
! ----- Compute strains
        epsi = 0.d0
        epsi = matmul(B, disp)
! ----- Compute stresses
        sigm(1+(kpg-1)*nbsig:nbsig*kpg) = matmul(hookeMatrix, epsi)
    enddo
!
end subroutine
