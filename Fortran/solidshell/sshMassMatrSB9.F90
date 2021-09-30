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
subroutine sshMassMatrSB9(nb_node , nb_dof  ,&
                          phenom  , j_mater , volume,&
                          matrMass)
!
implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/rcvalb.h"
#include "asterfort/elrefe_info.h"
!
integer, intent(in) :: nb_node, nb_dof
character(len=16), intent(in) :: phenom
integer, intent(in) :: j_mater
real(kind=8), intent(in) :: volume
real(kind=8), intent(out) :: matrMass(nb_dof, nb_dof)
!
! --------------------------------------------------------------------------------------------------
!
! Solid-shell element - SB9
!
! Compute mass matrix
!
! --------------------------------------------------------------------------------------------------
!
! In  nb_node          : number of nodes of element without pinch node(s)
! In  nb_dof           : number of dof
! In  phenom           : phenomenon to get material properties
! In  j_mater          : JEVEUX address for material properties
! In  volume           : volume of element
! Out matrMass         : mass matrix
!
! --------------------------------------------------------------------------------------------------
!
    integer :: jv_poids, jv_coopg, npg
    integer :: icodre(1)
    real(kind=8) ::  XI(3), poids
    integer :: i, j, L1, L2, C1, C2, kpg
    real(kind=8) :: g1(8), g2(8), g3(8)
    real(kind=8) :: h1(8), h2(8), h3(8), h4(8), s1(8)
    real(kind=8) :: N(8), N9
    real(kind=8) :: matrMassPt(25, 25)
    real(kind=8) :: rho(1)
    character(len=4), parameter :: fami = 'MASS'
    real(kind=8), parameter :: eye(3,3) = reshape((/1.d0, 0.d0, 0.d0,&
                                                    0.d0, 1.d0, 0.d0,&
                                                    0.d0, 0.d0, 1.d0/),(/3,3/))
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(nb_node .eq. 8)
    ASSERT(nb_dof .eq. 25)
    matrMass = 0.d0
!
! - Finite element informations
!
    call elrefe_info(fami=fami, npg=npg,&
                     jpoids=jv_poids, jcoopg=jv_coopg)
!
! - Vectors for decomposition
!
    s1 = (/ 1.d0,  1.d0,  1.d0,  1.d0,  1.d0,  1.d0, 1.d0,  1.d0/)/8.d0
    g1 = (/-1.d0,  1.d0,  1.d0, -1.d0, -1.d0,  1.d0, 1.d0, -1.d0/)/8.d0
    g2 = (/-1.d0, -1.d0,  1.d0,  1.d0, -1.d0, -1.d0, 1.d0,  1.d0/)/8.d0
    g3 = (/-1.d0, -1.d0, -1.d0, -1.d0,  1.d0,  1.d0, 1.d0,  1.d0/)/8.d0
    h1 = (/ 1.d0, -1.d0,  1.d0, -1.d0,  1.d0, -1.d0, 1.d0, -1.d0/)/8.d0
    h2 = (/ 1.d0,  1.d0, -1.d0, -1.d0, -1.d0, -1.d0, 1.d0,  1.d0/)/8.d0
    h3 = (/ 1.d0, -1.d0, -1.d0,  1.d0, -1.d0,  1.d0, 1.d0, -1.d0/)/8.d0
    h4 = (/-1.d0,  1.d0, -1.d0,  1.d0,  1.d0, -1.d0, 1.d0, -1.d0/)/8.d0
!
! - Loop on Gauss points
!
    do kpg = 1, npg
! ----- Current Gauss point
        poids = zr(jv_poids+kpg-1)
        XI(1) = zr(jv_coopg+3*(kpg-1)-1+1) 
        XI(2) = zr(jv_coopg+3*(kpg-1)-1+2)
        XI(3) = zr(jv_coopg+3*(kpg-1)-1+3)
! ----- Get density
        call rcvalb(fami, kpg   , 1  , '+'      , j_mater,&
                    ' ' , phenom, 0  , ' '      , [0.d0] ,&
                    1   , 'RHO' , rho, icodre(1), 1)
! ----- Construction of shape functions
        N  = s1+&
             XI(1)*g1+XI(2)*g2+XI(3)*g3+&
             XI(1)*XI(2)*h1+XI(3)*XI(2)*h2+XI(1)*XI(3)*h3+XI(1)*XI(2)*XI(3)*h4
        N9 = 1.d0 - XI(3)*XI(3)
! ----- Compute matrix on volumic nodes
        matrMassPt = 0.d0
        do i = 1, nb_node
            do j = 1, nb_node
                L1 = 3*(i-1)+1
                L2 = 3*(i-1)+3
                C1 = 3*(j-1)+1
                C2 = 3*(j-1)+3
                matrMassPt(L1:L2,C1:C2) = rho(1)*volume*poids*N(i)*N(j)*eye
            end do
        end do
! ----- Compute matrix on pinch node
        do i = 1, nb_node
            matrMassPt(nb_dof,3*(i-1)+3) = rho(1)*volume*poids*N(i)*N9
            matrMassPt(3*(i-1)+3,nb_dof) = rho(1)*volume*poids*N(i)*N9
        enddo
        matrMassPt(nb_dof, nb_dof) = rho(1)*poids*N9*N9
! ----- Add to matrix
        do i = 1, nb_dof
            do j = 1, nb_dof
                matrMass(i,j) = matrMass(i,j) + matrMassPt(i,j)
            end do
        end do 
    end do
!
end subroutine
