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
subroutine sshNLPreLogSB9(lVect   , lgpg    , vim ,&
                          epsgPrev, epsgCurr, &
                          epslPrev, epslIncr,&
                          gn      , lamb    , logl,&
                          tPrev   , iret)
!
implicit none
!
#include "asterf_types.h"
#include "asterfort/sshEpslSB9.h"
#include "blas/dcopy.h"
!
aster_logical, intent(in) :: lVect
integer, intent(in) :: lgpg
real(kind=8), intent(in) :: vim(lgpg)
real(kind=8), intent(in) :: epsgCurr(6), epsgPrev(6)
real(kind=8), intent(out) :: epslPrev(6), epslIncr(6)
real(kind=8), intent(out) :: gn(3, 3), lamb(3), logl(3)
real(kind=8), intent(out) :: tPrev(6)
integer, intent(out) :: iret
!
! --------------------------------------------------------------------------------------------------
!
! Solid-shell element
!
! Pre-treatment for GDEF_LOG
!
! --------------------------------------------------------------------------------------------------
!
! In  lVect            : flag to compute vector
! In  lgpg             : length of internal state variable vector
! In  vim              : internal state variables at beginning of time step
! In  epsgPrev         : Green-Lagrange strains at beginning of time step
! In  epsgCurr         : Green-Lagrange strains at end of time step
! Out gn               : eigen vectors for F tensor
! Out lamb             : eigen values for F tensor
! Out logl             : log(lamb)
! Out epslPrev         : logarithmic strains at beginning of time step
! Out epslIncr         : increment of logarithmic strains
! Out tPrev            : "logarithmic" stresses at beginning of time step
! Out iret             : return code for error
!                         0=OK, 1=vp(Ft.F) trop petites (compression infinie)
!
! --------------------------------------------------------------------------------------------------
!
    integer :: i, ivtn
    real(kind=8) :: epslCurr(6) 
!
! --------------------------------------------------------------------------------------------------
!
    gn       = 0.d0
    lamb     = 0.d0
    logl     = 0.d0
    epslPrev = 0.d0
    epslIncr = 0.d0
    tPrev    = 0.d0
    iret     = 0
!
! - Compute previous logarithmic strains
!
    call sshEpslSB9(epsgPrev, epslPrev, gn, lamb, logl, iret)
!
    if (iret .ne. 0) then
        goto 999
    endif
!
! - Compute incremental logarithmic strains
!
    if (lVect) then
! ----- Compute current logarithmic strains
        call sshEpslSB9(epsgCurr, epslCurr, gn, lamb, logl, iret)
        if (iret .ne. 0) then
            goto 999
        endif
! ----- Compute incremental logarithmic strains
        do i = 1, 6
            epslIncr(i) = epslCurr(i) - epslPrev(i)
        end do
    endif
!
! - Get "logarithmic" stresses from internal state variables at previous time step
!
    ivtn = lgpg-6+1
    call dcopy(6, vim(ivtn), 1, tPrev, 1)
!
999 continue
!
end subroutine
