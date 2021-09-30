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
subroutine sshNLStabSigmSB9(lLarge     ,&
                            nb_node    , nb_dof    ,&
                            geomInit   , disp      , Ueff,&
                            BXIdev     , BETAdev   ,&
                            BETAZETAdev, BXIZETAdev,&
                            SHXI       , SHETA     ,&
                            SHETAZETA  , SHXIZETA)
!
implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/sshGreenCartSB9.h"
!
aster_logical, intent(in) :: lLarge
integer, intent(in) :: nb_node, nb_dof
real(kind=8), intent(in)  :: geomInit(3*nb_node), disp(nb_dof), Ueff
real(kind=8), intent(in)  :: BXIdev(6,3*nb_node), BETAdev(6,3*nb_node)
real(kind=8), intent(in)  :: BXIZETAdev(6,3*nb_node), BETAZETAdev(6,3*nb_node)
real(kind=8), intent(out) :: SHXI(6), SHETA(6), SHETAZETA(6), SHXIZETA(6)
!
! --------------------------------------------------------------------------------------------------
!
! Solid-shell element - SB9
!
! Compute stabilisation stresses
!
! --------------------------------------------------------------------------------------------------
!
! In  lLarge           : flag for large strains
! In  nb_node          : number of nodes of element without pinch node(s)
! In  nb_dof           : number of dof
! In  geomInit         : initial coordinates of element
! In  disp             : displacements of element
! In  Ueff             : effective shear modulus
! In  BXIdev           : xi part of gradient matrix (deviatoric part)
! In  BETAdev          : eta part of gradient matrix (deviatoric part)
! In  BETAZETAdev      : eta x zeta part of gradient matrix (deviatoric part)
! In  BXIZETAdev       : xi x zeta part of gradient matrix (deviatoric part)
! Out SHXI             : stress stabilisation - xi x xi part
! Out SHETA            : stress stabilisation - eta x eta part
! Out SHETAZETA        : stress stabilisation - eta x zeta part
! Out SHXIZETA         : stress stabilisation - xi x zeta part
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: EH1(6), EH2(6), EH23(6), EH13(6), Dstab(6,6)
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(nb_node .eq. 8)
    ASSERT(nb_dof .eq. 25)
    SHXI      = 0.d0
    SHETA     = 0.d0
    SHETAZETA = 0.d0
    SHXIZETA  = 0.d0
!
! - Stabilisation matrix
!
    Dstab      =  0.d0
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
    Dstab = Dstab*Ueff
    !WRITE(6,*) 'Cont-Stab - DSTAB: ',sum(Dstab), Ueff
    !WRITE(6,*) 'Cont-Stab - Disp: ',sum(disp)
!
! - Compute the deformation gradient in cartesian base
!
    if (lLarge) then
        call sshGreenCartSB9(nb_node , nb_dof,&
                             geomInit, disp  ,&
                             EH1     , EH2   ,&
                             EH23    , EH13)
    else
        EH1  = matmul(BXIdev     , disp(1:3*nb_node))
        EH2  = matmul(BETAdev    , disp(1:3*nb_node))
        EH23 = matmul(BETAZETAdev, disp(1:3*nb_node))
        EH13 = matmul(BXIZETAdev , disp(1:3*nb_node))
    endif
!
! - Compute stabilisation stresses
!
    SHXI      = matmul(Dstab, EH1)
    SHETA     = matmul(Dstab, EH2)
    SHXIZETA  = matmul(Dstab, EH13)
    SHETAZETA = matmul(Dstab, EH23)
!
    !WRITE(6,*) 'Cont-Stab: ',SHXI,SHETA,SHXIZETA,SHETAZETA
!
end subroutine
