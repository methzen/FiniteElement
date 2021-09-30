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
subroutine sshNLStabMatrSB9(nb_node   , Ueff       ,&
                            BXIdev    , BETAdev    ,&
                            BXIZETAdev, BETAZETAdev,&
                            Kmstab)
!
implicit none
!
#include "asterfort/assert.h"
#include "asterfort/sshBTSB.h"
#include "asterfort/sshNLStabSigmSB9.h"
!
integer, intent(in) :: nb_node
real(kind=8), intent(in)  :: Ueff
real(kind=8), intent(in)  :: BXIdev(6,3*nb_node), BETAdev(6,3*nb_node)
real(kind=8), intent(in)  :: BXIZETAdev(6,3*nb_node), BETAZETAdev(6,3*nb_node)
real(kind=8), intent(out) :: Kmstab(3*nb_node,3*nb_node)
!
! --------------------------------------------------------------------------------------------------
!
! Solid-shell element - SB9
!
! Compute stabilization for material rigidity
!
! --------------------------------------------------------------------------------------------------
!
! In  nb_node          : number of nodes of element without pinch node(s)
! In  Ueff             : effective shear modulus
! In  BXIdev           : xi part of gradient matrix (deviatoric part)
! In  BETAdev          : eta part of gradient matrix (deviatoric part)
! In  BXIZETAdev       : xi x zeta part of gradient matrix (deviatoric part)
! In  BETAZETAdev      : eta x zeta part of gradient matrix (deviatoric part)
! Out Kmstab           : stabilization matrix for material rigidity
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8), dimension(6,6) :: Dstab
    real(kind=8), dimension(24,24) :: SXI, SETA, SXIZETA, SETAZETA
    real(kind=8) :: aux1, aux2
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(nb_node .eq. 8)
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
    call sshBTSB(Dstab, 6, 3*nb_node, BXIdev     , SXI)
    call sshBTSB(Dstab, 6, 3*nb_node, BETAdev    , SETA)
    call sshBTSB(Dstab, 6, 3*nb_node, BETAZETAdev, SETAZETA)
    call sshBTSB(Dstab, 6, 3*nb_node, BXIZETAdev , SXIZETA)
!
! - Final matrix
!
    Kmstab = 0.d0
    aux1   = 8.d0/3.d0
    aux2   = 8.d0/9.d0
    Kmstab = (SXI+SETA)*aux1+(SXIZETA+SETAZETA)*aux2
!
    !WRITE(6,*) 'Coefficient de stabilisation: ',Ueff
    !WRITE(6,*) 'Stabilization matrix for material rigidity: ',Kmstab
!
end subroutine
