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
subroutine sshBTSigmSB7(npg      , nbsig   ,&
                        nb_node  , nb_dof  ,&
                        jv_coopg , jv_poids,&
                        geom     , sigm    ,&
                        btSigm)
!
implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/btsig.h"
#include "asterfort/sshGradMatrSB7.h"
!
integer, intent(in) :: npg, nbsig, nb_node, nb_dof
integer, intent(in) :: jv_coopg, jv_poids
real(kind=8), intent(in) :: geom(3*nb_node), sigm(nbsig*npg)
real(kind=8), intent(out) :: btSigm(nb_dof)
!
! --------------------------------------------------------------------------------------------------
!
! Solid-shell element - SB7
!
! Compute B^T . SIGMA
!
! --------------------------------------------------------------------------------------------------
!
! In  npg              : number of Gauss points
! In  nbsig            : number stress components
! In  nb_node          : number of nodes of element without pinch node(s)
! In  nb_dof           : number of dof
! In  jv_coopg         : JEVEUX adress to coordinates of Gauss points
! In  jv_poids         : JEVEUX adress for weight of Gauss points 
! In  geom             : coordinates of element
! In  sigm             : stress at Gauss points
! Out btSigm           : product B^T . sigma
!
! --------------------------------------------------------------------------------------------------
!
    integer :: kpg
    real(kind=8) :: B(6, 19), poids, zeta, det
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(nb_node .eq. 6)
    ASSERT(nb_dof .eq. 19)
    ASSERT(nbsig .eq. 6)
    btSigm(1:nb_dof) = 0.d0
    do kpg = 1, npg
        poids = zr(jv_poids+kpg-1)
        zeta  = zr(jv_coopg-1+3*kpg)
! ----- Combute B matrix
        call sshGradMatrSB7(nb_node, nb_dof, geom, zeta, B, jacobian_ = det)
! ----- Product BT . sigma CUMULATED (see btsig subroutine)
        call btsig(nb_dof, nbsig, det*poids, B, sigm(1+nbsig*(kpg-1)), btSigm)
    end do
!
end subroutine
