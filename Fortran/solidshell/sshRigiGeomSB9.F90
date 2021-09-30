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
subroutine sshRigiGeomSB9(nbsig   , nb_node, nb_dof, npg, &
                          geom    , volume , sigm  ,&
                          matrGeom)
!
implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/sshShapeDeri.h"
#include "asterfort/sshTMatrSB9.h"
#include "asterfort/sshTMatrDecoSB9.h"
#include "asterfort/sshRigiGeomPtSB9.h"
#include "asterfort/elrefe_info.h"
!
integer, intent(in) :: nb_node, nb_dof, nbsig, npg
real(kind=8), intent(in) :: geom(24), volume
real(kind=8), intent(in) :: sigm(npg*nbsig)
real(kind=8), intent(out) :: matrGeom(nb_dof, nb_dof)
!
! --------------------------------------------------------------------------------------------------
!
! Solid-shell element - SB9
!
! Compute geometric matrix
!
! --------------------------------------------------------------------------------------------------
!
! In  nbsig            : number of components o the stress tensor
! In  nb_node          : number of nodes of element without pinch node(s)
! In  nb_dof           : number of dof
! In  npg              : number of integration points
! In  geom             : coordinates of element
! In  volume           : volume of current element
! In  sigm             : stresses
! Out matrGeom         : geometry rigidity matrix at current Gauss point
!
! --------------------------------------------------------------------------------------------------
!
    character(len=4), parameter :: fami = 'RIGI'
    integer :: jv_coopg, jv_poids
    integer :: kpg
    real(kind=8) :: zeta, poids
    real(kind=8) :: TZETA(6,6), T(6,6)
    real(kind=8) :: matrGeomPt(25, 25), sigmKpg(6)
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(nb_node .eq. 8)
    ASSERT(nb_dof .eq. 25)
    ASSERT(nbsig .eq. 6)
    matrGeom = 0.d0
!
! - Finite element informations
!
    call elrefe_info(fami=fami,&
                     jpoids=jv_poids, jcoopg=jv_coopg)
!
! - Compute the decomposition of the inverse Jacobian matrix
!
    call sshTMatrDecoSB9(nb_node, geom, TZETA_ = TZETA)
!
! - Compute T matrix (matrix relating the covariant and cartesian frames) at center of element
!
    call sshTMatrSB9(nb_node, geom, (/0.d0,0.d0,0.d0/), T)
!
! - Loop on Gauss points
!
    do kpg = 1, npg
! ----- Current Gauss point
        poids = zr(jv_poids+kpg-1)
        zeta  = zr(jv_coopg-1+3*kpg)
! ----- Compute geometric matrix
        sigmKpg = sigm(1+(kpg-1)*nbsig:nbsig*kpg)
        call sshRigiGeomPtSB9(nbsig     , nb_node,&
                              nb_dof    , &
                              TZETA     , T      ,&
                              zeta      , sigmKpg,&
                              matrGeomPt)
!
        matrGeom = matrGeom + poids*volume*matrGeomPt
    end do
!
end subroutine
