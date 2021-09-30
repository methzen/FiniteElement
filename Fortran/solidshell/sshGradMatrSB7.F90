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
subroutine sshGradMatrSB7(nb_node, nb_dof, geom     , zeta,&
                          B      , BStab_, jacobian_)
!
implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/matinv.h"
#include "asterfort/sshGradMatrCovaSB7.h"
#include "asterfort/sshGradMatrPinchSB7.h"
#include "asterfort/sshGradMatrMembSB7.h"
#include "asterfort/sshJaco.h"
#include "asterfort/sshLocalFrameSB7.h"
#include "asterfort/sshGradMatrBendSB7.h"
#include "asterfort/sshGradMatrTwistSB7.h"
#include "asterfort/sshGradMatrShearSB7.h"
!
integer, intent(in) :: nb_node, nb_dof
real(kind=8), intent(in) :: geom(3*nb_node), zeta
real(kind=8), intent(out) :: B(6,19)
real(kind=8), optional, intent(out) :: BStab_(2,18), jacobian_
!
! --------------------------------------------------------------------------------------------------
!
! Solid-shell element - SB7
!
! Compute B matrix at current Gauss point
!
! --------------------------------------------------------------------------------------------------
!
! In  nb_node          : number of nodes of element without pinch node(s)
! In  nb_dof           : number of dof
! In  zeta             : out-of-plane parametric component
! In  geom             : coordinates of element
! Out B                : gradient B matrix
! Out BStab            : pinch part stabilization of the gradient B matrix
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: geomLocal(18), R(3,3)
    real(kind=8) :: h, triaArea, z
    real(kind=8) :: bx(6), by(6), bz(6)
    real(kind=8) :: triaN1x, triaN1y, triaN2x, triaN2y, triaN3x, triaN3y
    real(kind=8) :: Bm0(3,18), Bb0(3,18), Bp0(18)
    real(kind=8) :: Bc0(2,18), Bc1(2,18), Bc2(2,18)
    real(kind=8) :: J1s3(3,3), invJ1s3(3,3)
    real(kind=8), parameter :: center(3)  = (/1.d0/3.d0, 1.d0/3.d0, 0.d0/)
    real(kind=8) :: BStab(2,18), jacobian
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(nb_node .eq. 6)
    ASSERT(nb_dof .eq. 19)
    B        = 0.d0
    BStab    = 0.d0
    jacobian = 0.d0
!
! - Compute local frame axis, local coordinates and vertices of mid_edge triangle
!
    call sshLocalFrameSB7(geom    ,&
                          R       , geomLocal, h,&
                          triaN1x , triaN1y  ,&
                          triaN2x , triaN2y  ,&
                          triaN3x , triaN3y  ,&
                          triaArea)
    !WRITE(6,*) 'h: ',h
    !WRITE(6,*) 'geom: ',sum(geom)
    !WRITE(6,*) 'geomLocal: ',sum(geomLocal)
    !WRITE(6,*) 'R: ',sum(R)
!
! - Compute gradient matrix in covariant basis at center of prism
!
    call sshGradMatrCovaSB7(nb_node, geomLocal, center, bx, by, bz)
    !WRITE(6,*) 'Bx: ',sum(bx),sum(by),sum(bz)
!
! - Compute membrane terms of the gradient B matrix
!
    call sshGradMatrMembSB7(R, bx, by, Bm0)
    !WRITE(6,*) 'Bm0: ',sum(Bm0)
!
! - Compute pinching terms of the gradient B matrix
!
    call sshGradMatrPinchSB7(R, bz, Bp0)
    !WRITE(6,*) 'Bp0: ',sum(Bp0)
!
! - Compute bending terms of the gradient B matrix
!
    call sshGradMatrBendSB7(nb_node , geom   , R,&
                            triaN1x , triaN1y,&
                            triaN2x , triaN2y,&
                            triaN3x , triaN3y,&
                            triaArea, Bb0)
    !WRITE(6,*) 'Bb0: ',sum(Bb0)
!
! - Compute shear terms of the gradient B matrix
!
    call sshGradMatrShearSB7(geomLocal, R, Bc0, Bc1, Bc2)
    !WRITE(6,*) 'Bc0: ',sum(Bc0)
    !WRITE(6,*) 'Bc1: ',sum(Bc1)
    !WRITE(6,*) 'Bc2: ',sum(Bc2)
!
! - Compute B matrix at current Gauss point
!
    z = zeta*h/2.d0
    B(1:2,1:18) = Bm0(1:2,1:18)+z*Bb0(1:2,1:18)
    B(3,  1:18) = Bp0
    B(4,  1:18) = Bm0(3,1:18)+z*Bb0(3,1:18)
    B(5:6,1:18) = (Bc0+Bc1/3.d0+Bc2/3.d0)*(1.d0-zeta*zeta)*1.25d0
    B(  3,  19) = -4.d0*zeta/h
    !WRITE(6,*) 'B: ',sum(B), jacobian
!
! - Compute term for stabilization of twist mode
!
    call sshGradMatrTwistSB7(geomLocal, R, BStab)
    !WRITE(6,*) 'BStab: ',sum(BStab)
!
! - Compute jacobian
!
    call sshJaco(nb_node, geom, (/1.d0/3.d0, 1.d0/3.d0, zeta/), J1s3)
    call matinv('S', 3, J1s3, invJ1s3, jacobian)
!
    if (present(jacobian_)) then
        jacobian_ = jacobian
    endif
!
    if (present(BStab_)) then
        BStab_ = BStab
    endif
!
end subroutine
