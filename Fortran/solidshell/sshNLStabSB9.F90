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
! aslint: disable=W1306
!
subroutine sshNLStabSB9(lMatr     , lVect      , lMatrPred, lLarge,&
                        nb_node   , nb_dof     ,&
                        geomInit  , dispPrev   , dispCurr,&
                        det       , Ueff       ,&
                        BXIdev    , BETAdev    ,&
                        BXIZETAdev, BETAZETAdev,&
                        matuu     , vectu)
!
implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/sshNLStabMatrSB9.h"
#include "asterfort/sshNLStabVectSB9.h"
#include "asterfort/sshNLStabMatrGeomSB9.h"
!
aster_logical, intent(in) :: lMatr, lVect, lMatrPred, lLarge
integer, intent(in) :: nb_node, nb_dof
real(kind=8), intent(in) :: geomInit(3*nb_node), dispCurr(nb_dof), dispPrev(nb_dof)
real(kind=8), intent(in)  :: Ueff, det
real(kind=8), intent(in) :: BXIdev(6,24), BETAdev(6,24), BETAZETAdev(6,24), BXIZETAdev(6,24)
real(kind=8), intent(inout) :: vectu(*), matuu(*)
!
! --------------------------------------------------------------------------------------------------
!
! Solid-shell element
!
! Compute stabilization for matrix and internal force
!
! --------------------------------------------------------------------------------------------------
!
! In  lMatr            : flag to compute matrix
! In  lVect            : flag to compute vector
! In  lMatrPred        : flag to compute matrix when prediction
! In  lLarge           : flag for large strains
! In  nb_node          : number of nodes of element without pinch node(s)
! In  nb_dof           : number of dof
! In  geomInit         : initial coordinates of element
! In  dispPrev         : displacements of element at beginning of step time
! In  dispCurr         : displacements of element at end of step time
! In  det              : determinant of transformation
! In  Ueff             : effective shear modulus
! In  BXIdev           : xi part of gradient matrix 
! In  BETAdev          : eta part of gradient matrix
! In  BXIZETAdev       : xi x zeta part of gradient matrix
! In  BETAZETAdev      : eta x zeta part of gradient matrix
! IO  vectu            : linear internal force vector
! IO  matuu            : elastic matrix
!
! --------------------------------------------------------------------------------------------------
!
    integer :: i, j, k
    real(kind=8) :: Kgstab(24,24), Kmstab(3*nb_node,3*nb_node), Fstab(3*nb_node)
    real(kind=8) :: disp(nb_dof)
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(nb_node .eq. 8)
    ASSERT(nb_dof .eq. 25)
!
! - Select displacement
!
    if (lMatrPred) then
        disp = dispPrev
    else
        disp = dispCurr
    endif
!
! - Matrix
!
    if (lMatr) then
! ----- Material part
        call sshNLStabMatrSB9(nb_node   , Ueff       ,&
                              BXIdev    , BETAdev    ,&
                              BXIZETAdev, BETAZETAdev,&
                              Kmstab)
        k = 0
        do i = 1, 3*nb_node
            do j = 1, i
                k = k + 1
                matuu(k) = matuu(k) + det*Kmstab(i,j)
            enddo
        enddo
        !WRITE(6,*) 'Matrice de stabilisation matérielle: ',sum(Kmstab)
! ----- Geometric part
        if (lLarge) then
            call sshNLStabMatrGeomSB9(nb_node    , nb_dof,&
                                      geomInit   , disp      , Ueff,&
                                      BXIdev     , BETAdev   ,&
                                      BETAZETAdev, BXIZETAdev,&
                                      Kgstab)
            k = 0
            do i = 1, 3*nb_node
                do j = 1, i
                    k = k + 1
                    matuu(k) = matuu(k) + det*Kgstab(i,j)
                enddo
            enddo
            !WRITE(6,*) 'Matrice de stabilisation géométrique: ',sum(Kgstab)
        endif
    endif 
!
! - Vector
!
    if (lVect) then
        call sshNLStabVectSB9(lLarge     ,&
                              nb_node    , nb_dof    ,&
                              geomInit   , dispCurr  , Ueff,&
                              BXIdev     , BETAdev   ,&
                              BETAZETAdev, BXIZETAdev,&
                              Fstab)
        do i = 1, 3*nb_node
            vectu(i) = vectu(i) + det*Fstab(i)
        enddo
    endif 
!
end subroutine
