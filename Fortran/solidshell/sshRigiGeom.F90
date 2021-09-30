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
subroutine sshRigiGeom()
!
implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/jevech.h"
#include "asterfort/assert.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/matinv.h"
#include "asterfort/nbsigm.h"
#include "asterfort/sshJaco.h"
#include "asterfort/sshRigiGeomSB7.h"
#include "asterfort/sshRigiGeomSB9.h"
!
! --------------------------------------------------------------------------------------------------
!
! Solid-shell element
!
! Compute geometric matrix - Option RIGI_GEOM
!
! --------------------------------------------------------------------------------------------------
!
    integer, parameter :: nbDofMaxi = 25, nbNodeMaxi = 8
    integer :: i, j, k
    integer :: jv_matr, jv_geom, jv_sigm
    integer :: nb_dof, nb_node, npg, nno, nbsig
    real(kind=8) :: matrGeom(nbDofMaxi, nbDofMaxi)
    real(kind=8) :: geomInit(3*nbNodeMaxi), center(3)
    real(kind=8) :: J0(3,3), Jm1(3,3), volume
    character(len=4), parameter :: fami = 'RIGI'
    aster_logical :: sb9, sb7
!
! --------------------------------------------------------------------------------------------------
!
    nbsig = nbsigm()
    ASSERT(nbsig .eq. 6)
    matrGeom = 0.d0
!
! - Finite element informations
!
    call elrefe_info(fami=fami, nno=nno, npg=npg)
!
! - Sizes
!
    nb_node = 0
    sb7     = ASTER_FALSE
    sb9     = ASTER_FALSE
    if (nno .eq. 9) then
        nb_node = 8
        sb9     = ASTER_TRUE
        center  = (/0.d0, 0.d0, 0.d0/)
    elseif (nno .eq. 7) then
        nb_node = 6
        sb7     = ASTER_TRUE
        center  = (/1.d0/3.d0, 1.d0/3.d0, 0.d0/)
    else
        ASSERT(ASTER_FALSE)
    endif
    nb_dof = 3*(nb_node) + 1
    ASSERT(nb_dof .le. nbDofMaxi)
    ASSERT(nb_node .le. nbNodeMaxi)
!
! - Geometry
!
    call jevech('PGEOMER', 'L', jv_geom)
!
! - Get initial geometry
!
    do i = 1, 3*nb_node
        geomInit(i) = zr(jv_geom+i-1)
    end do
!
! - Stresses
!
    call jevech('PCONTRR', 'L', jv_sigm)
!
! - Compute initial jacobian matrix at center of element
!
    call sshJaco(nb_node, geomInit, center, J0)
!
! - Inverse jacobian matrix
!
    call matinv('S', 3, J0, Jm1, volume)
!
! - Misoriented cell => get absolute !
!
    volume=abs(volume)
!
! - Compute geometric rigidity
!
    if (sb9) then
        call sshRigiGeomSB9(nbsig   , nb_node, nb_dof     , npg, &
                            geomInit, volume , zr(jv_sigm),&
                            matrGeom)
    elseif (sb7) then
        call sshRigiGeomSB7(nbsig   , nb_node, nb_dof     , npg, &
                            geomInit, volume , zr(jv_sigm),&
                            matrGeom)
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
            zr(jv_matr+k-1) = matrGeom(i,j)
        end do
    end do
!
end subroutine
