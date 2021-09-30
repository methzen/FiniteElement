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
subroutine sshTwistVectSB7(geomLocal, invJo, Vgamma)
!
implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/matinv.h"
!
real(kind=8), intent(in) :: geomLocal(18)
real(kind=8), intent(out) :: Vgamma(6, 2), invJo(3, 3)
!
! --------------------------------------------------------------------------------------------------
!
! Solid-shell element - SB7
!
! Compute twist vector for stabilization
!
! --------------------------------------------------------------------------------------------------
!
! In  geomLocal        : coordinates of element in the frame of mid_edge triangle
! Out invJo            : inverse of jacobian in ksi=0, eta=0, zeta=0
! Out Vgamma           : twist vector for stabilization
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8), parameter :: uns2 = 1.d0/2.d0
    integer :: i, j
    real(kind=8) :: XiL, YiL, ZiL, XjL, YjL, ZjL, XkL, YkL, ZkL
    real(kind=8) :: XlL, YlL, ZlL, XmL, YmL, ZmL, XnL, YnL, ZnL
    real(kind=8) :: hX, hY, hZ
    real(kind=8) :: detJo, J0(3,3)
    real(kind=8) :: Bksi(3,6), Bx(6), By(6), Bz(6)
    real(kind=8) :: XYZe(6,3), matH(6,2)
!
! --------------------------------------------------------------------------------------------------
!
    Vgamma = 0.d0
    invJo  = 0.d0
!
! - Get local coordinates of nodes of prism in the frame of mid_edge triangle
!
    XiL = geomLocal(1 )
    YiL = geomLocal(2 )
    ZiL = geomLocal(3 )
    XjL = geomLocal(4 )
    YjL = geomLocal(5 )
    ZjL = geomLocal(6 )
    XkL = geomLocal(7 )
    YkL = geomLocal(8 )
    ZkL = geomLocal(9 )
    XlL = geomLocal(10)
    YlL = geomLocal(11)
    ZlL = geomLocal(12)
    XmL = geomLocal(13)
    YmL = geomLocal(14)
    ZmL = geomLocal(15)
    XnL = geomLocal(16)
    YnL = geomLocal(17)
    ZnL = geomLocal(18)
!
! - Compute jacobien in ksi=0, eta=0, zeta=0
!
    J0(1,1)=(-XiL+XjL   -XlL+XmL   )/2.d0
    J0(1,2)=(-YiL+YjL   -YlL+YmL   )/2.d0
    J0(1,3)=(-ZiL+ZjL   -ZlL+ZmL   )/2.d0
    J0(2,1)=(-XiL   +XkL-XlL   +XnL)/2.d0
    J0(2,2)=(-YiL   +YkL-YlL   +YnL)/2.d0
    J0(2,3)=(-ZiL   +ZkL-ZlL   +ZnL)/2.d0
    J0(3,1)=(-XiL      +XlL        )/2.d0
    J0(3,2)=(-YiL      +YlL        )/2.d0
    J0(3,3)=(-ZiL      +ZlL        )/2.d0
    call matinv('S', 3, J0, invJo, detJo)
!
! - Compute Bksi matrix
!
    Bksi(1,1)=-uns2
    Bksi(2,1)=-uns2
    Bksi(3,1)=-uns2
    Bksi(1,2)= uns2
    Bksi(2,2)= 0.d0
    Bksi(3,2)= 0.d0
    Bksi(1,3)= 0.d0
    Bksi(2,3)= uns2
    Bksi(3,3)= 0.d0
    Bksi(1,4)=-uns2
    Bksi(2,4)=-uns2
    Bksi(3,4)= uns2
    Bksi(1,5)= uns2
    Bksi(2,5)= 0.d0
    Bksi(3,5)= 0.d0
    Bksi(1,6)= 0.d0
    Bksi(2,6)= uns2
    Bksi(3,6)= 0.d0
!
    do j = 1, 6
        Bx(j) = invJo(1,1)*Bksi(1,j)+invJo(1,2)*Bksi(2,j)+invJo(1,3)*Bksi(3,j)
        By(j) = invJo(2,1)*Bksi(1,j)+invJo(2,2)*Bksi(2,j)+invJo(2,3)*Bksi(3,j)
        Bz(j) = invJo(3,1)*Bksi(1,j)+invJo(3,2)*Bksi(2,j)+invJo(3,3)*Bksi(3,j)
    enddo
!
! - [Vgamma]=VecteursGamma(xyz,Bx,By,Bz)
!
    XYZe(1,1)=XiL
    XYZe(1,2)=YiL
    XYZe(1,3)=ZiL
    XYZe(2,1)=XjL
    XYZe(2,2)=YjL
    XYZe(2,3)=ZjL
    XYZe(3,1)=XkL
    XYZe(3,2)=YkL
    XYZe(3,3)=ZkL
    XYZe(4,1)=XlL
    XYZe(4,2)=YlL
    XYZe(4,3)=ZlL
    XYZe(5,1)=XmL
    XYZe(5,2)=YmL
    XYZe(5,3)=ZmL
    XYZe(6,1)=XnL
    XYZe(6,2)=YnL
    XYZe(6,3)=ZnL
!
! - vecteurs {h}
!
    matH(1,1)= 0.d0
    matH(1,2)= 0.d0
    matH(2,1)= 0.d0
    matH(2,2)=-1.d0
    matH(3,1)=-1.d0
    matH(3,2)= 0.d0
    matH(4,1)= 0.d0
    matH(4,2)= 0.d0
    matH(5,1)= 0.d0
    matH(5,2)= 1.d0
    matH(6,1)= 1.d0
    matH(6,2)= 0.d0
!
! - Compute Vgamma
!
    do j = 1, 2
        hX = 0.d0
        hY = 0.d0
        hZ = 0.d0
        do i = 1, 6
            hX = hX+matH(i,j)*XYZe(i,1)
            hY = hY+matH(i,j)*XYZe(i,2)
            hZ = hZ+matH(i,j)*XYZe(i,3)
        enddo
        do i = 1, 6
            Vgamma(i,j) = (matH(i,j)-hX*Bx(i)-hY*By(i)-hZ*Bz(i))/2.d0
        enddo
    enddo
!
end subroutine
