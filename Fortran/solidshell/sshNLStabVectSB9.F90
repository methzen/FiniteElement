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
subroutine sshNLStabVectSB9(lLarge     ,&
                            nb_node    , nb_dof    ,&
                            geomInit   , disp      , Ueff,&
                            BXIdev     , BETAdev   ,&
                            BETAZETAdev, BXIZETAdev,&
                            Fstab)
!
implicit none
!
#include "asterfort/assert.h"
#include "asterfort/sshNLStabSigmSB9.h"
!
aster_logical, intent(in) :: lLarge
integer, intent(in) :: nb_node, nb_dof
real(kind=8), intent(in)  :: geomInit(3*nb_node), disp(nb_dof), Ueff
real(kind=8), intent(in)  :: BXIdev(6,3*nb_node), BETAdev(6,3*nb_node)
real(kind=8), intent(in)  :: BXIZETAdev(6,3*nb_node), BETAZETAdev(6,3*nb_node)
real(kind=8), intent(out) :: Fstab(3*nb_node)
!
! --------------------------------------------------------------------------------------------------
!
! Solid-shell element - SB9
!
! Compute stabilisation for vector of internal forces
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
! Out Fstab            : stabilization for internal forces
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: aux1, aux2
    real(kind=8) :: SHXI(6), SHETA(6), SHETAZETA(6), SHXIZETA(6)
!
! --------------------------------------------------------------------------------------------------
!
    Fstab = 0.d0
!
! - Compute stabilisation stresses
!
    call sshNLStabSigmSB9(lLarge     ,&
                          nb_node    , nb_dof    ,&
                          geomInit   , disp      , Ueff,&
                          BXIdev     , BETAdev   ,&
                          BETAZETAdev, BXIZETAdev,&
                          SHXI       , SHETA     ,&
                          SHETAZETA  , SHXIZETA)
!
! - Compute stabilisation for internal forces
!
    aux1  = 8.d0/3.d0
    aux2  = 8.d0/9.d0
    Fstab = (matmul(transpose(BXIdev),SHXI) +&
             matmul(transpose(BETAdev),SHETA))*aux1+&
            (matmul(transpose(BETAZETAdev),SHETAZETA)+&
             matmul(transpose(BXIZETAdev),SHXIZETA))*aux2
!
end subroutine
