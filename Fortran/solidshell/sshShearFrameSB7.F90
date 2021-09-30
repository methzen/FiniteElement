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
! aslint: disable=W1504
!
subroutine sshShearFrameSB7(XiL , YiL , ZiL , XjL, YjL, ZjL,&
                            XkL , YkL , ZkL , XlL, YlL, ZlL,&
                            XmL , YmL , ZmL , XnL, YnL, ZnL,&
                            g304, g305, g306,&
                            V1  , V2  , V3  ,&
                            h1  , h2  , h3)
!
implicit none
!
#include "asterfort/assert.h"
!
real(kind=8), intent(in) :: XiL, YiL, ZiL, XjL, YjL, ZjL, XkL, YkL, ZkL
real(kind=8), intent(in) :: XlL, YlL, ZlL, XmL, YmL, ZmL, XnL, YnL, ZnL
real(kind=8), intent(out) :: g304(3), g305(3), g306(3)
real(kind=8), intent(out) :: V1(3,3), V2(3,3), V3(3,3)
real(kind=8), intent(out) :: h1, h2, h3
!
! --------------------------------------------------------------------------------------------------
!
! Solid-shell element - SB7
!
! Compute frames of SB7 element to evaluate transversal shear
!
! --------------------------------------------------------------------------------------------------
!
! In  XiL, YiL, ZiL    : local coordinates in the frame of mid_edge triangle for node 1 of prism
! In  XjL, YjL, ZjL    : local coordinates in the frame of mid_edge triangle for node 2 of prism
! In  XkL, YkL, ZkL    : local coordinates in the frame of mid_edge triangle for node 3 of prism
! In  XlL, YlL, ZlL    : local coordinates in the frame of mid_edge triangle for node 4 of prism
! In  XmL, YmL, ZmL    : local coordinates in the frame of mid_edge triangle for node 5 of prism
! In  XnL, YnL, ZnL    : local coordinates in the frame of mid_edge triangle for node 6 of prism
! Out g304             : natural frame for first edge of mid-edge triangle
! Out g305             : natural frame for second edge of mid-edge triangle
! Out g306             : natural frame for third edge of mid-edge triangle
! Out V1               : orthogonal frame defined at first node of mid-edge triangle
! Out V2               : orthogonal frame defined at second node of mid-edge triangle
! Out V3               : orthogonal frame defined at third node of mid-edge triangle
! Out h1               : length of "vertical" edge at first node of mid-edge triangle
! Out h2               : length of "vertical" edge at second node of mid-edge triangle
! Out h3               : length of "vertical" edge at third node of mid-edge triangle
!
! --------------------------------------------------------------------------------------------------
!
    integer :: i
    real(kind=8) :: X10(3), X20(3), X30(3)
    real(kind=8) :: Nx, Ny
    real(kind=8) :: V3xx, V1xx, V2xx, V2yx, V2zx, V3yx, V3zx, V1yx, V1zx
    real(kind=8) :: V3xy, V1xy, V2xy, V2yy, V2zy, V3yy, V3zy, V1yy, V1zy
    real(kind=8) :: V3xz, V1xz, V2xz, V2yz, V2zz, V3yz, V3zz, V1yz, V1zz
    real(kind=8) :: x1, x2, x3, y1, y2, y3, z1, z2, z3
!
! --------------------------------------------------------------------------------------------------
!

!
! - Normalized length of edge at vertexes of mid-edge triangle
!
    h1     = sqrt((XlL-XiL)*(XlL-XiL)+(YlL-YiL)*(YlL-YiL)+(ZlL-ZiL)*(ZlL-ZiL))
    X10(1) = (XlL-XiL)/h1
    X10(2) = (YlL-YiL)/h1
    X10(3) = (ZlL-ZiL)/h1
    h2     = sqrt((XmL-XjL)*(XmL-XjL)+(YmL-YjL)*(YmL-YjL)+(ZmL-ZjL)*(ZmL-ZjL))
    X20(1) = (XmL-XjL)/h2
    X20(2) = (YmL-YjL)/h2
    X20(3) = (ZmL-ZjL)/h2
    h3     = sqrt((XnL-XkL)*(XnL-XkL)+(YnL-YkL)*(YnL-YkL)+(ZnL-ZkL)*(ZnL-ZkL))
    X30(1) = (XnL-XkL)/h3
    X30(2) = (YnL-YkL)/h3
    X30(3) = (ZnL-ZkL)/h3
!
! - Natural frame for each edge of mid-edge triangle
!
    do i = 1, 3
        g304(i) = (h1*X10(i)+h2*X20(i))/4.d0
        g305(i) = (h2*X20(i)+h3*X30(i))/4.d0
        g306(i) = (h3*X30(i)+h1*X10(i))/4.d0
    enddo
!
! - Coordinates of point for middle of each edge of mid-edge triangle
!
    x1 = (XiL+XlL)/2.d0
    y1 = (YiL+YlL)/2.d0
    z1 = (ZiL+ZlL)/2.d0
    x2 = (XjL+XmL)/2.d0
    y2 = (YjL+YmL)/2.d0
    z2 = (ZjL+ZmL)/2.d0
    x3 = (XkL+XnL)/2.d0
    y3 = (YkL+YnL)/2.d0
    z3 = (ZkL+ZnL)/2.d0
!
! - Orthogonal frame defined at first node of mid-edge triangle
!
    V1zx = (XlL-XiL)/h1
    V1zy = (YlL-YiL)/h1
    V1zz = (ZlL-ZiL)/h1
    V1yx = V1zy*(z2-z1)-V1zz*(y2-y1)
    V1yy = V1zz*(x2-x1)-V1zx*(z2-z1)
    V1yz = V1zx*(y2-y1)-V1zy*(x2-x1)
    Ny   = sqrt(V1yx*V1yx+V1yy*V1yy+V1yz*V1yz)
    V1yx = V1yx/Ny
    V1yy = V1yy/Ny
    V1yz = V1yz/Ny
    V1xx = V1yy*V1zz-V1zy*V1yz
    V1xy = V1yz*V1zx-V1zz*V1yx
    V1xz = V1yx*V1zy-V1zx*V1yy
    Nx   = sqrt(V1xx*V1xx+V1xy*V1xy+V1xz*V1xz)
    V1xx = V1xx/Nx
    V1xy = V1xy/Nx
    V1xz = V1xz/Nx
    V1(1,1:3) = (/V1xx,V1xy,V1xz/)
    V1(2,1:3) = (/V1yx,V1yy,V1yz/)
    V1(3,1:3) = (/V1zx,V1zy,V1zz/)
!
! - Orthogonal frame defined at second node of mid-edge triangle
!
    V2zx = (XmL-XjL)/h2
    V2zy = (YmL-YjL)/h2
    V2zz = (ZmL-ZjL)/h2
    V2yx = V2zy*(z3-z2)-V2zz*(y3-y2)
    V2yy = V2zz*(x3-x2)-V2zx*(z3-z2)
    V2yz = V2zx*(y3-y2)-V2zy*(x3-x2)
    Ny   = sqrt(V2yx*V2yx+V2yy*V2yy+V2yz*V2yz)
    V2yx = V2yx/Ny
    V2yy = V2yy/Ny
    V2yz = V2yz/Ny
    V2xx = V2yy*V2zz-V2zy*V2yz
    V2xy = V2yz*V2zx-V2zz*V2yx
    V2xz = V2yx*V2zy-V2zx*V2yy
    Nx   = sqrt(V2xx*V2xx+V2xy*V2xy+V2xz*V2xz)
    V2xx = V2xx/Nx
    V2xy = V2xy/Nx
    V2xz = V2xz/Nx
    V2(1,1:3) = (/V2xx,V2xy,V2xz/)
    V2(2,1:3) = (/V2yx,V2yy,V2yz/)
    V2(3,1:3) = (/V2zx,V2zy,V2zz/)
!
! - Orthogonal frame defined at third node of mid-edge triangle
!
    V3zx = (XnL-XkL)/h3
    V3zy = (YnL-YkL)/h3
    V3zz = (ZnL-ZkL)/h3
    V3yx = V3zy*(z1-z3)-V3zz*(y1-y3)
    V3yy = V3zz*(x1-x3)-V3zx*(z1-z3)
    V3yz = V3zx*(y1-y3)-V3zy*(x1-x3)
    Ny   = sqrt(V3yx*V3yx+V3yy*V3yy+V3yz*V3yz)
    V3yx = V3yx/Ny
    V3yy = V3yy/Ny
    V3yz = V3yz/Ny
    V3xx = V3yy*V3zz-V3zy*V3yz
    V3xy = V3yz*V3zx-V3zz*V3yx
    V3xz = V3yx*V3zy-V3zx*V3yy
    Nx   = sqrt(V3xx*V3xx+V3xy*V3xy+V3xz*V3xz)
    V3xx = V3xx/Nx
    V3xy = V3xy/Nx
    V3xz = V3xz/Nx
    V3(1,1:3) = (/V3xx,V3xy,V3xz/)
    V3(2,1:3) = (/V3yx,V3yy,V3yz/)
    V3(3,1:3) = (/V3zx,V3zy,V3zz/)
!
end subroutine
