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
! but WITHOUT ANY WARRANTY without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with code_aster.  If not, see <http://www.gnu.org/licenses/>.
! --------------------------------------------------------------------
!
subroutine sshStabMatrSB9(nb_node, Ueff    ,&
                          BXI    , BETA    ,&
                          BXIZETA, BETAZETA,&
                          SXI    , SETA    ,&
                          SXIZETA, SETAZETA)
!
implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/sshBTSB.h"
!
integer, intent(in) :: nb_node
real(kind=8), intent(in) :: Ueff
real(kind=8), intent(in) :: BXI(6,3*nb_node), BETA(6,3*nb_node)
real(kind=8), intent(in) :: BXIZETA(6,3*nb_node), BETAZETA(6,3*nb_node)
real(kind=8), intent(out) :: SXI(3*nb_node,3*nb_node), SETA(3*nb_node,3*nb_node)
real(kind=8), intent(out) :: SETAZETA(3*nb_node,3*nb_node), SXIZETA(3*nb_node,3*nb_node)
!
! --------------------------------------------------------------------------------------------------
!
! Solid-shell element - SB9
!
! Compute stabilization matrix (elastic)
!
! --------------------------------------------------------------------------------------------------
!
! In  nb_node          : number of nodes of element without pinch node(s)
! In  Ueff             : elastic shear modulus
! In  BXI              : xi part of gradient matrix
! In  BETA             : eta part of gradient matrix
! In  BXIZETA          : xi x zeta part of gradient matrix
! In  BETAZETA         : eta x zeta part of gradient matrix
! Out SXI              : stabilization matrix - xi x xi part
! Out SETA             : stabilization matrix - eta x eta part
! Out SXIZETA          : stabilization matrix - xi x zeta part
! Out SETAZETA         : stabilization matrix - eta x zeta part
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: Dstab(6,6)
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(nb_node .eq. 8)
!
! - Compute the deformation gradient B in cartesian base - Deviatoric part
! - D'apres Hughes ces termes sont deja deviatoric dans leur forme initiale.
!

!
! - Compute stabilization
!
    Dstab      = 0.d0
    Dstab(1,1) = +4.d0/3.d0
    Dstab(1,2) = -2.d0/3.d0
    Dstab(1,3) = -2.d0/3.d0
    Dstab(2,1) = -2.d0/3.d0
    Dstab(2,2) = +4.d0/3.d0
    Dstab(2,3) = -2.d0/3.d0
    Dstab(3,1) = -2.d0/3.d0
    Dstab(3,2) = -2.d0/3.d0
    Dstab(3,3) = +4.d0/3.d0
    Dstab(4,4) = 1.d0
    Dstab(5,5) = 1.d0
    Dstab(6,6) = 1.d0
    Dstab      = Dstab*Ueff
!
! - Products
!
    call sshBTSB(Dstab, 6, 3*nb_node, BXI     , SXI)
    call sshBTSB(Dstab, 6, 3*nb_node, BETA    , SETA)
    call sshBTSB(Dstab, 6, 3*nb_node, BETAZETA, SETAZETA)
    call sshBTSB(Dstab, 6, 3*nb_node, BXIZETA , SXIZETA)
!
end subroutine
