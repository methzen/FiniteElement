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
subroutine sshRigiGeomPtSB7(nbsig     , nb_node  , nb_dof,&
                            h         , geomLocal, R,&
                            zeta      , sigm     ,&
                            matrGeomPt)
!
implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/sshBTSB.h"
#include "asterfort/sshGradMatrCovaSB7.h"
!
integer, intent(in) :: nb_node, nb_dof, nbsig
real(kind=8), intent(in) :: zeta, h, sigm(nbsig)
real(kind=8), intent(in)  :: R(3, 3), geomLocal(18)
real(kind=8), intent(out) :: matrGeomPt(nb_dof, nb_dof)
!
! --------------------------------------------------------------------------------------------------
!
! Solid-shell element - SB7
!
! Compute geometric matrix at Gauss point
!
! --------------------------------------------------------------------------------------------------
!
! In  nbsig            : number of components o the stress tensor
! In  nb_node          : number of nodes of element without pinch node(s)
! In  nb_dof           : number of dof
! In  geomLocal        : coordinates of element in the frame of mid_edge triangle
! In  R                : rotation matrix for basis of mid_edge triangle
! In  h                : element thickness
! In  zeta             : zeta (out-of-plane parametric coordinate) of current Gauss point
! In  sigm             : stresses of current Gauss point
! Out matrGeomPt       : geometry rigidity matrix at current Gauss point
!
! --------------------------------------------------------------------------------------------------
!
    integer :: i, j
    real(kind=8) :: XI(3)
    real(kind=8) :: bx(6), by(6), bz(6)
    real(kind=8) :: BG(9, 18), BGTot(9, 19)
    real(kind=8) :: SHEAD(9,9), SIG(3,3)
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(nb_node .eq. 6)
    ASSERT(nb_dof .eq. 19)
    ASSERT(nbsig .eq. 6)
    matrGeomPt = 0.d0
    !WRITE(6,*) 'MatrGeom-Sigm: ', sum(sigm)
    !WRITE(6,*) 'MatrGeom-Rot: ',sum(R)
    !WRITE(6,*) 'MatrGeom-GeomLocal: ',sum(geomLocal), h
!
! - Compute gradient matrix in covariant basis at center of prism
!
    XI = (/1.d0/3.d0, 1.d0/3.d0, zeta/)
    call sshGradMatrCovaSB7(nb_node, geomLocal, XI, bx, by, bz)
    !WRITE(6, *) 'MatrGeom-XI: ',XI
    !WRITE(6, *) 'MatrGeom-Bx: ',sum(bx),sum(by), sum(bz)
!
! - Compute B matrix for geometry
!
    BG = 0.d0
    do j = 1, 6
          BG(1,3*(j-1)+1  ) = R(1,1)*bx(j)
          BG(1,3*(j-1)+2  ) = R(1,2)*bx(j)
          BG(1,3*(j-1)+3  ) = R(1,3)*bx(j)
          BG(2,3*(j-1)+1  ) = R(1,1)*by(j)
          BG(2,3*(j-1)+2  ) = R(1,2)*by(j)
          BG(2,3*(j-1)+3  ) = R(1,3)*by(j)
          BG(3,3*(j-1)+1  ) = R(1,1)*bz(j)
          BG(3,3*(j-1)+2  ) = R(1,2)*bz(j)
          BG(3,3*(j-1)+3  ) = R(1,3)*bz(j)
          BG(3+1,3*(j-1)+1) = R(2,1)*bx(j)
          BG(3+1,3*(j-1)+2) = R(2,2)*bx(j)
          BG(3+1,3*(j-1)+3) = R(2,3)*bx(j)
          BG(3+2,3*(j-1)+1) = R(2,1)*by(j)
          BG(3+2,3*(j-1)+2) = R(2,2)*by(j)
          BG(3+2,3*(j-1)+3) = R(2,3)*by(j)
          BG(3+3,3*(j-1)+1) = R(2,1)*bz(j)
          BG(3+3,3*(j-1)+2) = R(2,2)*bz(j)
          BG(3+3,3*(j-1)+3) = R(2,3)*bz(j)
          BG(6+1,3*(j-1)+1) = R(3,1)*bx(j)
          BG(6+1,3*(j-1)+2) = R(3,2)*bx(j)
          BG(6+1,3*(j-1)+3) = R(3,3)*bx(j)
          BG(6+2,3*(j-1)+1) = R(3,1)*by(j)
          BG(6+2,3*(j-1)+2) = R(3,2)*by(j)
          BG(6+2,3*(j-1)+3) = R(3,3)*by(j)
          BG(6+3,3*(j-1)+1) = R(3,1)*bz(j)
          BG(6+3,3*(j-1)+2) = R(3,2)*bz(j)
          BG(6+3,3*(j-1)+3) = R(3,3)*bz(j)
    enddo
    !WRITE(6, *) 'Geom-BG sans pinch: ',sum(BG)
!
! - Taking additional node in consideration
!
    BGtot = 0.d0
    do I = 1, 9
        do J = 1, 18
            BGtot(i,j) = BG(i,j)
        end do
    end do
    BGtot(9,19) = -4.d0*zeta/h
    !WRITE(6, *) 'Geom-BG avec pinch: ',sum(BGtot)
!
! - Matrix for stress
!
    SIG(1,1:3) = (/sigm(1),0.d0,0.d0/)
    SIG(2,1:3) = (/0.d0,sigm(2),0.d0/)
    SIG(3,1:3) = (/0.d0,0.d0,sigm(3)/)
    !WRITE(6,*) 'Geom-SIG: ',sum(SIG)
    SHEAD = 0.d0
    do I = 1, 3
        do J = 1, 3
            SHEAD(i,j)     = SIG(i,j)
            SHEAD(i+3,j+3) = SIG(i,j)
            SHEAD(i+6,j+6) = SIG(i,j)
        end do
    end do
    !WRITE(6,*) 'Geom-SHEAD: ',sum(SHEAD)
!
! - Compute geometric rigidity
!
    call sshBTSB(SHEAD, 9, 19, BGtot, matrGeomPt)
    !WRITE(6,*) 'Geom-matrGeomPt: ',sum(matrGeomPt)
!
end subroutine
