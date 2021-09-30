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
subroutine sshNLPosLogSB9(lVect , lMatr  , lgpg,&
                          tPrev , tCurr  ,&
                          gn    , lamb   , logl,&
                          dtde  , vip    , &
                          dsidep, pk2Curr)
!
implicit none
!
#include "asterf_types.h"
#include "asterfort/deflg2.h"
#include "asterfort/deflg3.h"
#include "asterfort/pmat.h"
#include "asterfort/symt46.h"
#include "blas/daxpy.h"
#include "blas/dcopy.h"
!
aster_logical, intent(in) :: lVect, lMatr
integer, intent(in) :: lgpg
real(kind=8), intent(in) :: tPrev(6), tCurr(6)
real(kind=8), intent(in) :: gn(3, 3), lamb(3), logl(3)
real(kind=8), intent(in) :: dtde(6, 6)
real(kind=8), intent(inout) :: vip(lgpg)
real(kind=8), intent(out) :: dsidep(6, 6), pk2Curr(6)
!
! --------------------------------------------------------------------------------------------------
!
! Solid-shell element - SB9
!
! Post-treatment for GDEF_LOG
!
! --------------------------------------------------------------------------------------------------
!
! In  lVect            : flag to compute vector
! In  lMatr            : flag to compute matrix
! In  lgpg             : length of internal state variable vector
! In  tPrev            : "logarithmic" stresses at beginning of time step
! In  tCurr            : "logarithmic" stresses at end of time step
! In  gn               : eigen vectors for F tensor
! In  lamb             : eigen values for F tensor
! In  logl             : log(lamb)
! In  dtde             : tangent matrix (dT / dE)
! IO  vip              : internal state variables at end of time step
! Out dsidep           : tangent matrix after transformation
! Out pk2Curr          : stresses (PK2) at end of time step
!
! --------------------------------------------------------------------------------------------------
!
    integer :: i, j, ivtn
    real(kind=8) :: pes_t(6, 6), trav2(6, 6)
    real(kind=8) :: pes(6, 6),  tWork(6)
    real(kind=8) :: tl(3, 3, 3, 3), tl_syme(6, 6)
    real(kind=8) :: feta(4), xi(3, 3), me(3, 3, 3, 3)
!
! --------------------------------------------------------------------------------------------------
!
    dsidep  = 0.d0
    pk2Curr = 0.d0
!
! - Compute tensor P (symmetric)
!
    call deflg2(gn, lamb, logl, pes, feta, xi, me)
    pes_t = transpose(pes)
!
! - Compute tangent matrix
!
    if (lMatr) then
! ----- Get "logarithmic" stresses T
        if (lVect) then
            call dcopy(6, tCurr, 1, tWork, 1)
        else
            call dcopy(6, tPrev, 1, tWork, 1)
        endif
! ----- Compute tensor T:L
        call deflg3(gn, feta, xi, me, tWork, tl)
! ----- Symmetric version of T:L
        call symt46(tl, tl_syme)
! ----- Compute dsidep
        call pmat(6, pes_t, dtde, trav2)
        call pmat(6, trav2, pes, dsidep)
        call daxpy(36, 1.d0, tl_syme, 1, dsidep,1)
    endif
!
! - Compute PK2 stresses
!
    if (lVect) then
        do i = 1, 6
            do j = 1, 6
                pk2Curr(i) = pk2Curr(i) + tCurr(j)*pes(j,i)
            end do
        end do
    endif
!
! - Save "logarithmic" stresses in internal state variables
!
    if (lVect) then
        ivtn = lgpg-6+1
        call dcopy(6, tCurr, 1, vip(ivtn), 1)
    endif
end subroutine
