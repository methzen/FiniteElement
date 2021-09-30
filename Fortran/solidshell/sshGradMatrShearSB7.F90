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
subroutine sshGradMatrShearSB7(geomLocal, R, Bc0, Bc1, Bc2)
!
implicit none
!
real(kind=8), intent(in)  :: geomLocal(18), R(3,3)
real(kind=8), intent(out) :: Bc1(2,18), Bc2(2,18), Bc0(2,18)
!
#include "asterfort/matinv.h"
#include "asterfort/sshShearFrameSB7.h"
!
! --------------------------------------------------------------------------------------------------
!
! Solid-shell element - SB7
!
! Compute shear terms of the gradient B matrix
!
! --------------------------------------------------------------------------------------------------
!
! In  geomLocal        : coordinates of element in the frame of mid_edge triangle
! In  R                : rotation matrix for frame of mid_edge triangle
! Out Bcx              : shear parts of the gradient B matrix
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: XiL, YiL, ZiL, XjL, YjL, ZjL, XkL, YkL, ZkL
    real(kind=8) :: XlL, YlL, ZlL, XmL, YmL, ZmL, XnL, YnL, ZnL
    real(kind=8) :: g304(3), g305(3), g306(3)
    real(kind=8) :: V1(3,3), V2(3,3), V3(3,3)
    real(kind=8) :: h1, h2, h3
    real(kind=8) :: x1,y1,x2,y2,x3,y3
    real(kind=8) :: triaNormX,triaNormY,triaNormZ
    real(kind=8) :: triaTau2X,triaTau2Y,triaTau2Z
    real(kind=8) :: triaTau1X,triaTau1Y,triaTau1Z
    integer      :: q,p,rr,i
    real(kind=8) :: g10(3),g20(3),g21(3)
    real(kind=8) :: C111(5),C112(5),C212(5),C213(5),C311(5),C313(5),C113(5)
    real(kind=8) :: C211(5),C312(5),C1(2,15),C2(2,15),C3(2,15)
    real(kind=8) :: C121(5),D1(2,2),D2(2,2),D3(2,2)
    real(kind=8) :: C122(5)
    real(kind=8) :: C123(5)
    real(kind=8) :: C221(5)
    real(kind=8) :: C222(5)
    real(kind=8) :: C223(5),Bdt1(2,18)
    real(kind=8) :: C321(5),Bdt2(2,18)
    real(kind=8) :: C322(5),Bdt3(2,18)
    real(kind=8) :: C323(5),z1,z2,z3
    real(kind=8) :: Bd1(2,15),aux1,T(15,18)
    real(kind=8) :: Bd2(2,15),aux2,Tm(6,18)
    real(kind=8) :: Bd3(2,15),aux3,Tr(6,18)
    real(kind=8) :: detJo, J0(3,3), invJo(3,3)
!
! --------------------------------------------------------------------------------------------------
!
    Bc0 = 0.d0
    Bc1 = 0.d0
    Bc2 = 0.d0
!
! - Get basis for mid-edge triangle
!
    triaNormX = R(1,1)
    triaNormY = R(1,2)
    triaNormZ = R(1,3)
    triaTau1X = R(3,1)
    triaTau1Y = R(3,2)
    triaTau1Z = R(3,3)
    triaTau2X = R(2,1)
    triaTau2Y = R(2,2)
    triaTau2Z = R(2,3)
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
! - Coordinates of nodes for mid-edge triangle
!
    x1 = (XiL+XlL)/2.d0
    y1 = (YiL+YlL)/2.d0
    x2 = (XjL+XmL)/2.d0
    y2 = (YjL+YmL)/2.d0
    x3 = (XkL+XnL)/2.d0
    y3 = (YkL+YnL)/2.d0
    z1 = (ZiL+ZlL)/2.d0
    z2 = (ZjL+ZmL)/2.d0
    z3 = (ZkL+ZnL)/2.d0
!
! - Coordinates of nodes for mid-edge triangle with node 1 as origin
!
    x3 = x3-x1
    y3 = y3-y1
    x2 = x2-x1
    y2 = 0.d0
    x1 = 0.d0
    y1 = 0.d0
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
! - Compute frames of SB7 element to evaluate transversal shear
!
    call sshShearFrameSB7(XiL , YiL , ZiL , XjL, YjL, ZjL,&
                          XkL , YkL , ZkL , XlL, YlL, ZlL,&
                          XmL , YmL , ZmL , XnL, YnL, ZnL,&
                          g304, g305, g306,&
                          V1  , V2  , V3  ,&
                          h1  , h2  , h3)
!
! - Vectors for edges of mid-edge triangle
!
    g10(1) = x2-x1
    g10(2) = y2-y1
    g10(3) = z2-z1
    g20(1) = x3-x1
    g20(2) = y3-y1
    g20(3) = z3-z1
    g21(1) = x3-x2
    g21(2) = y3-y2
    g21(3) = z3-z2
!
! - Compute C
!
    C111 = 0.d0
    C112 = 0.d0
    C212 = 0.d0
    C213 = 0.d0
    C311 = 0.d0
    C313 = 0.d0
    C113 = 0.d0
    C211 = 0.d0
    C312 = 0.d0
    C121 = 0.d0
    C222 = 0.d0
    C122 = 0.d0
    C123 = 0.d0
    C221 = 0.d0
    C223 = 0.d0
    C321 = 0.d0
    C322 = 0.d0
    C323 = 0.d0
    C1   = 0.d0
    C2   = 0.d0
    C3   = 0.d0
    C111(1)=-g304(1)
    C111(2)=-g304(2)
    C111(3)=-g304(3)
    C111(4)=-0.25d0*h1*sum(V1(2,1:3)*g10)
    C111(5)= 0.25d0*h1*sum(V1(1,1:3)*g10)
!
    C112(1)=+g304(1)
    C112(2)=+g304(2)
    C112(3)=+g304(3)
    C112(4)=-0.25d0*h2*sum(V2(2,1:3)*g10)
    C112(5)=+0.25d0*h2*sum(V2(1,1:3)*g10)
!
    C212(1)=-g305(1)
    C212(2)=-g305(2)
    C212(3)=-g305(3)
    C212(4)=-0.25d0*h2*sum(V2(2,1:3)*g21)
    C212(5)=+0.25d0*h2*sum(V2(1,1:3)*g21)
!
    C213(1)=+g305(1)
    C213(2)=+g305(2)
    C213(3)=+g305(3)
    C213(4)=-0.25d0*h3*sum(V3(2,1:3)*g21)
    C213(5)=+0.25d0*h3*sum(V3(1,1:3)*g21)
!
    C311(1)=+g306(1)
    C311(2)=+g306(2)
    C311(3)=+g306(3)
    C311(4)=+0.25d0*h1*sum(V1(2,1:3)*g20)
    C311(5)=-0.25d0*h1*sum(V1(1,1:3)*g20)
!
    C313(1)=-g306(1)
    C313(2)=-g306(2)
    C313(3)=-g306(3)
    C313(4)=+0.25d0*h3*sum(V3(2,1:3)*g20)
    C313(5)=-0.25d0*h3*sum(V3(1,1:3)*g20)
!
    do i = 1, 5
        C121(i)=-C311(i)
        C122(i)=-C312(i)
        C123(i)=-C313(i)
        C221(i)=-C111(i)
        C222(i)=-C112(i)
        C223(i)=-C113(i)
        C321(i)=-C211(i)
        C322(i)=-C212(i)
        C323(i)=-C213(i)
    end do
    do i = 1, 3
        C1(1,i)  =C111(i)
        C1(1,i+3)=C112(i)
        C1(1,i+6)=C113(i)
        C1(2,i)  =C121(i)
        C1(2,i+3)=C122(i)
        C1(2,i+6)=C123(i)
        C2(1,i)  =C211(i)
        C2(1,i+3)=C212(i)
        C2(1,i+6)=C213(i)
        C2(2,i)  =C221(i)
        C2(2,i+3)=C222(i)
        C2(2,i+6)=C223(i)
        C3(1,i)  =C311(i)
        C3(1,i+3)=C312(i)
        C3(1,i+6)=C313(i)
        C3(2,i)  =C321(i)
        C3(2,i+3)=C322(i)
        C3(2,i+6)=C323(i)
    enddo
!
    do i = 4, 5
        C1(1,i+6) =C111(i)
        C1(1,i+8) =C112(i)
        C1(1,i+10)=C113(i)
        C1(2,i+6) =C121(i)
        C1(2,i+8) =C122(i)
        C1(2,i+10)=C123(i)
        C2(1,i+6) =C211(i)
        C2(1,i+8) =C212(i)
        C2(1,i+10)=C213(i)
        C2(2,i+6) =C221(i)
        C2(2,i+8) =C222(i)
        C2(2,i+10)=C223(i)
        C3(1,i+6) =C311(i)
        C3(1,i+8) =C312(i)
        C3(1,i+10)=C313(i)
        C3(2,i+6) =C321(i)
        C3(2,i+8) =C322(i)
        C3(2,i+10)=C323(i)
    enddo
!
! - Compute D
!
    D1   = 0.d0
    D2   = 0.d0
    D3   = 0.d0
    D1(1,1)= 1.d0
    D1(1,2)= 0.d0
    D1(2,1)= 0.d0
    D1(2,2)= 1.d0
    D2(1,1)= 0.d0
    D2(1,2)=-1.d0
    D2(2,1)= 1.d0
    D2(2,2)=-1.d0
    D3(1,1)=-1.d0
    D3(1,2)= 1.d0
    D3(2,1)=-1.d0
    D3(2,2)= 0.d0
!
! - Compute Bd = D*C
!
    Bd1=0.d0
    Bd2=0.d0
    Bd3=0.d0
    do p = 1,2
        do q = 1,15
            aux1=0.d0
            aux2=0.d0
            aux3=0.d0
            do rr = 1,2
                aux1=aux1+D1(p,rr)*C1(rr,q)
                aux2=aux2+D2(p,rr)*C2(rr,q)
                aux3=aux3+D3(p,rr)*C3(rr,q)
            enddo
            Bd1(p,q)=aux1
            Bd2(p,q)=aux2
            Bd3(p,q)=aux3
        enddo
    enddo
!
! - Compute matrix T
!
    T = 0.d0
    do p = 1, 3
        T(3*p-2,3*p-2)=0.5d0*triaNormX
        T(3*p-2,3*p-1)=0.5d0*triaNormY
        T(3*p-2,3*p)  =0.5d0*triaNormZ
        T(3*p-1,3*p-2)=0.5d0*triaTau2X
        T(3*p-1,3*p-1)=0.5d0*triaTau2Y
        T(3*p-1,3*p)  =0.5d0*triaTau2Z
        T(3*p,3*p-2)  =0.5d0*triaTau1X
        T(3*p,3*p-1)  =0.5d0*triaTau1Y
        T(3*p,3*p)    =0.5d0*triaTau1Z
        T(3*p-2,3*p+7)=0.5d0*triaNormX
        T(3*p-2,3*p+8)=0.5d0*triaNormY
        T(3*p-2,3*p+9)=0.5d0*triaNormZ
        T(3*p-1,3*p+7)=0.5d0*triaTau2X
        T(3*p-1,3*p+8)=0.5d0*triaTau2Y
        T(3*p-1,3*p+9)=0.5d0*triaTau2Z
        T(3*p,3*p+7)  =0.5d0*triaTau1X
        T(3*p,3*p+8)  =0.5d0*triaTau1Y
        T(3*p,3*p+9)  =0.5d0*triaTau1Z
    enddo
!
    Tm=0.d0
    Tr=0.d0
    do i = 1,3
        Tm(1,i)  =+V1(2,i)/h1
        Tm(2,i)  =-V1(1,i)/h1
        Tm(1,i+9)=-Tm(1,i)
        Tm(2,i+9)=-Tm(2,i)

        Tm(3,i+3) =+V2(2,i)/h2
        Tm(4,i+3) =-V2(1,i)/h2
        Tm(3,i+12)=-Tm(3,i+3)
        Tm(4,i+12)=-Tm(4,i+3)

        Tm(5,i+6) =+V3(2,i)/h3
        Tm(6,i+6) =-V3(1,i)/h3
        Tm(5,i+15)=-Tm(5,i+6)
        Tm(6,i+15)=-Tm(6,i+6)
    enddo
!
    do i = 1, 2
        Tr(i,1)=Tm(i,1)*triaNormX+Tm(i,2)*triaTau2X+Tm(i,3)*triaTau1X
        Tr(i,2)=Tm(i,1)*triaNormY+Tm(i,2)*triaTau2Y+Tm(i,3)*triaTau1Y
        Tr(i,3)=Tm(i,1)*triaNormZ+Tm(i,2)*triaTau2Z+Tm(i,3)*triaTau1Z
        Tr(i,10)=-Tr(i,1)
        Tr(i,11)=-Tr(i,2)
        Tr(i,12)=-Tr(i,3)
    enddo
    do i = 3, 4
        Tr(i,4)=Tm(i,4)*triaNormX+Tm(i,5)*triaTau2X+Tm(i,6)*triaTau1X
        Tr(i,5)=Tm(i,4)*triaNormY+Tm(i,5)*triaTau2Y+Tm(i,6)*triaTau1Y
        Tr(i,6)=Tm(i,4)*triaNormZ+Tm(i,5)*triaTau2Z+Tm(i,6)*triaTau1Z
        Tr(i,13)=-Tr(i,4)
        Tr(i,14)=-Tr(i,5)
        Tr(i,15)=-Tr(i,6)
    enddo
    do i = 5, 6
        Tr(i,7)=Tm(i,7)*triaNormX+Tm(i,8)*triaTau2X+Tm(i,9)*triaTau1X
        Tr(i,8)=Tm(i,7)*triaNormY+Tm(i,8)*triaTau2Y+Tm(i,9)*triaTau1Y
        Tr(i,9)=Tm(i,7)*triaNormZ+Tm(i,8)*triaTau2Z+Tm(i,9)*triaTau1Z
        Tr(i,16)=-Tr(i,7)
        Tr(i,17)=-Tr(i,8)
        Tr(i,18)=-Tr(i,9)
    enddo
!
    do p = 1, 6
        do q = 1, 18
            T(p+9,q)=Tr(p,q)
        enddo
    enddo
!
    Bdt1=0.d0
    Bdt2=0.d0
    Bdt3=0.d0
    do p = 1, 2
        do q = 1, 18
            aux1=0.d0
            aux2=0.d0
            aux3=0.d0
            do rr = 1, 15
                aux1=aux1+Bd1(p,rr)*T(rr,q)
                aux2=aux2+Bd2(p,rr)*T(rr,q)
                aux3=aux3+Bd3(p,rr)*T(rr,q)
            enddo
            Bdt1(p,q)=aux1
            Bdt2(p,q)=aux2
            Bdt3(p,q)=aux3
        enddo
    enddo
!
    do p = 1, 2
        do q = 1, 18
            aux1=0.d0
            aux2=0.d0
            aux3=0.d0
            do rr = 1, 2
                aux1=aux1+invJo(p,rr)*Bdt1(rr,q)
                aux2=aux2+invJo(p,rr)*(Bdt3(rr,q)-Bdt1(rr,q))
                aux3=aux3+invJo(p,rr)*(Bdt2(rr,q)-Bdt1(rr,q))
            enddo
            Bc0(p,q)=aux1
            Bc1(p,q)=aux2
            Bc2(p,q)=aux3
        enddo
    enddo
!
end subroutine
