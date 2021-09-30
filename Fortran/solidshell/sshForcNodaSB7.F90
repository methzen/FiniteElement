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
subroutine sshForcNodaSB7(npg     , nbsig   ,&
                          nb_node , nb_dof  ,&
                          jv_poids, jv_coopg,&
                          geom    , sigm    ,&
                          btSigm)
!
implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/sshBTSigmSB7.h"
!
integer, intent(in) :: npg, nb_node, nb_dof, nbsig
integer, intent(in) :: jv_coopg, jv_poids
real(kind=8), intent(in) :: geom(3*nb_node), sigm(nbsig*npg)
real(kind=8), intent(out) :: btSigm(nb_dof)
!
! --------------------------------------------------------------------------------------------------
!
! Solid-shell element - SB7
!
! Compute FORC_NODA
!
! --------------------------------------------------------------------------------------------------
!
! In  npg              : number of Gauss points
! In  nbsig            : number of components o the stress tensor
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
    ASSERT(nb_node .eq. 6)
    ASSERT(nb_dof .eq. 19)
    ASSERT(nbsig .eq. 6)
    btSigm = 0.d0
!
! - Compute B^T . Sig
!
    call sshBTSigmSB7(npg      , nbsig   ,&
                      nb_node  , nb_dof  ,&
                      jv_coopg , jv_poids,&
                      geom     , sigm    ,&
                      btSigm)
!
end subroutine
