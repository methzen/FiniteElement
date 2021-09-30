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
subroutine sshMassMeca()
!
implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/jevech.h"
#include "asterfort/assert.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/matinv.h"
#include "asterfort/rccoma.h"
#include "asterfort/sshJaco.h"
#include "asterfort/sshMassMatrSB9.h"
!
! --------------------------------------------------------------------------------------------------
!
! Solid-shell element
!
! Compute mass matrix - Option Mass_MECA
!
! --------------------------------------------------------------------------------------------------
!
    integer :: icodre(1)
    integer :: j_mater, jv_mate, jv_matr, jv_geom
    integer :: nb_dof, nb_node, nno
    real(kind=8) :: matrMass(25, 25)
    character(len=16) :: phenom
    integer :: i, j, k
    real(kind=8) :: geomInit(24), center(3)
    real(kind=8) :: J0(3,3), Jm1(3,3), volume
    character(len=4), parameter :: fami = 'MASS'
    aster_logical :: sb9
!
! --------------------------------------------------------------------------------------------------
!
    matrMass = 0.d0
!
! - Finite element informations
!
    call elrefe_info(fami=fami, nno=nno)
!
! - Sizes
!
    nb_node = 0
    sb9     = ASTER_FALSE
    center  = 0.d0
    if (nno .eq. 9) then
        nb_node = 8
        sb9     = ASTER_TRUE
        center  = (/0.d0, 0.d0, 0.d0/)
    else
        ASSERT(ASTER_FALSE)
    endif
    nb_dof = 3*(nb_node) + 1
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
! - Material parameters
!
    call jevech('PMATERC', 'L', jv_mate)
    j_mater = zi(jv_mate)
!
! - Get properties
!
    call rccoma(j_mater, 'ELAS', 1, phenom, icodre(1))
!
! - Compute jacobian matrix at center of element
!
    call sshJaco(nb_node, geomInit, center , J0)
!
! - Inverse jacobian matrix
!
    call matinv('S', 3, J0, Jm1, volume)
!
! - Misoriented cell => get absolute !
!
    volume=abs(volume)
!
! - Compute mass matrix
!
    if (sb9) then
        call sshMassMatrSB9(nb_node , nb_dof ,&
                            phenom  , j_mater, volume,&
                            matrMass)
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
            zr(jv_matr+k-1) = matrMass(i,j)
        end do
    end do
!
end subroutine
