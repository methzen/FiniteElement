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
subroutine sshSiefElga()
!
implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterc/r8vide.h"
#include "asterfort/jevech.h"
#include "asterfort/assert.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/nbsigm.h"
#include "asterfort/sshSiefElgaSB7.h"
#include "asterfort/sshSiefElgaSB9.h"
!
! --------------------------------------------------------------------------------------------------
!
! Solid-shell element
!
! Compute stresses - Option SIEF_EGA
!
! --------------------------------------------------------------------------------------------------
!
    integer, parameter :: nbNpgMaxi = 5, nbDofMaxi = 25, nbNodeMaxi = 8
    integer :: i
    integer :: jv_poids, jv_coopg
    integer :: j_mater, jv_mate, jv_sigm, jv_geom, jv_disp
    integer :: nb_dof, nb_node, npg, nno, nbsig
    real(kind=8) :: sigm(6*nbNpgMaxi)
    real(kind=8) :: geomInit(24)
    character(len=4), parameter :: fami = 'RIGI'
    aster_logical :: sb9, sb7
!
! --------------------------------------------------------------------------------------------------
!
    nbsig = nbsigm()
    ASSERT(nbsig .eq. 6)
!
! - Finite element informations
!
    call elrefe_info(fami=fami, npg=npg, nno=nno,&
                     jpoids=jv_poids, jcoopg=jv_coopg)
    ASSERT(npg .le. nbNpgMaxi)
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
    ASSERT(nb_dof .le. nbDofMaxi)
    ASSERT(nb_node .le. nbNodeMaxi)
!
! - Geometry
!
    call jevech('PGEOMER', 'L', jv_geom)
    do i = 1, 3*nb_node
        geomInit(i) = zr(jv_geom+i-1)
    enddo
!
! - Material parameters
!
    call jevech('PMATERC', 'L', jv_mate)
    j_mater = zi(jv_mate)
!
! - Displacements
!
    call jevech('PDEPLAR', 'L', jv_disp)
!
! - Compute rigidity matrix
!
    if (sb9) then
        call sshSiefElgaSB9(npg     , nbsig      ,&
                            nb_node , nb_dof     ,&
                            j_mater , jv_coopg   ,&
                            geomInit, zr(jv_disp),&
                            sigm)
    elseif (sb7) then
        call sshSiefElgaSB7(npg     , nbsig      ,&
                            nb_node , nb_dof     ,&
                            j_mater , jv_coopg   ,&
                            geomInit, zr(jv_disp),&
                            sigm)
    else
        ASSERT(ASTER_FALSE)
    endif
!
! - Save stress
!
    call jevech('PCONTRR', 'E', jv_sigm)
    do i = 1, nbsig*npg
        zr(jv_sigm+i-1) = sigm(i)
    enddo
!
end subroutine
