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
subroutine sshRigiMeca()
!
implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterc/r8vide.h"
#include "asterfort/jevech.h"
#include "asterfort/assert.h"
#include "asterfort/sshRigiMatrSB7.h"
#include "asterfort/sshRigiMatrSB9.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/dmat3d.h"
!
! --------------------------------------------------------------------------------------------------
!
! Solid-shell element
!
! Compute rigidity matrix (linear) - Option RIGI_MECA
!
! --------------------------------------------------------------------------------------------------
!
    integer :: i, j, k
    integer :: jv_poids, jv_coopg
    integer :: j_mater, jv_mate, jv_matr, jv_geom
    integer :: nb_dof, nb_node, npg, nno
    real(kind=8) :: hookeMatrix(6, 6), matrRigi(25, 25)
    real(kind=8) :: repere(7), xyzgau(3)
    character(len=4), parameter :: fami = 'RIGI'
    aster_logical :: sb9, sb7
!
! --------------------------------------------------------------------------------------------------
!
    matrRigi    = 0.d0
    hookeMatrix = 0.d0
    repere      = 0.d0
    xyzgau      = 0.d0
!
! - Finite element informations
!
    call elrefe_info(fami=fami, npg=npg, nno=nno,&
                     jpoids=jv_poids, jcoopg=jv_coopg)
!
! - Sizes
!
    nb_node = 0
    sb7     = ASTER_FALSE
    sb9     = ASTER_FALSE
    if (nno .eq. 9) then
        nb_node = 8
        sb9     = ASTER_TRUE
    elseif (nno .eq. 7) then
        nb_node = 6
        sb7     = ASTER_TRUE
    else
        ASSERT(ASTER_FALSE)
    endif
    nb_dof = 3*(nb_node) + 1
!
! - Geometry
!
    call jevech('PGEOMER', 'L', jv_geom)
!
! - Material parameters
!
    call jevech('PMATERC', 'L', jv_mate)
    j_mater = zi(jv_mate)
!
! - Get Hooke matrix
!
    call dmat3d(fami, j_mater, r8vide(), '+'        , 1,&
                1   , repere , xyzgau  , hookeMatrix)
!
! - Compute rigidity matrix
!
    if (sb9) then
        call sshRigiMatrSB9(npg        , nb_node , nb_dof  ,&
                            jv_geom    , jv_coopg, jv_poids,&
                            hookeMatrix, matrRigi)
    elseif (sb7) then
        call sshRigiMatrSB7(npg        , nb_node , nb_dof  ,&
                            jv_geom    , jv_coopg, jv_poids,&
                            hookeMatrix, matrRigi)
    else
        ASSERT(ASTER_FALSE)
    endif
!
! - Save matrix
!
    call jevech('PMATUUR', 'E', jv_matr)
    k = 0
    do i = 1, nb_dof
        do j = 1, i
            k = k + 1
            zr(jv_matr+k-1) = matrRigi(i,j)
        end do
    end do
!
end subroutine
