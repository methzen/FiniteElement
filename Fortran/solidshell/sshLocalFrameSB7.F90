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
subroutine sshLocalFrameSB7(geom     ,&
                            R        , geomLocal, h,&
                            triaN1x_ , triaN1y_ ,&
                            triaN2x_ , triaN2y_ ,&
                            triaN3x_ , triaN3y_ ,&
                            triaArea_)
!
implicit none
!
#include "asterfort/assert.h"
!
real(kind=8), intent(in)  :: geom(18)
real(kind=8), intent(out) :: geomLocal(18), R(3,3), h
real(kind=8), optional, intent(out) :: triaN1x_, triaN1y_
real(kind=8), optional, intent(out) :: triaN2x_, triaN2y_
real(kind=8), optional, intent(out) :: triaN3x_, triaN3y_
real(kind=8), optional, intent(out) :: triaArea_
!
! --------------------------------------------------------------------------------------------------
!
! Solid-shell element - SB7
!
! Compute local frame axis, local coordinates and vertices of mid_edge triangle
!
! --------------------------------------------------------------------------------------------------
!
! The local frame is the frame of mid_edge triangle
!
! In  geom             : coordinates of element
! Out geomLocal        : coordinates of element in the frame of mid_edge triangle
! Out R                : rotation matrix for basis of mid_edge triangle
! Out h                : element thickness
! Out triaN1x, triaN1y : local coordinates of mid_edge triangle, first vertex
! Out triaN2x, triaN2y : local coordinates of mid_edge triangle, second vertex
! Out triaN3x, triaN3y : local coordinates of mid_edge triangle, third vertex
! Out triaArea         : area of mid_edge triangle
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: triaNormX, triaNormY, triaNormZ
    real(kind=8) :: Nx, Ny, Nz, X13, Y13, Z13, X12, Y12, Z12
    real(kind=8) :: triaTau2X, triaTau2Y, triaTau2Z
    real(kind=8) :: triaTau1X, triaTau1Y, triaTau1Z
    real(kind=8) :: triaArea
    real(kind=8) :: triaN1x, triaN1y, triaN2x, triaN2y, triaN3x, triaN3y
!
! --------------------------------------------------------------------------------------------------
!
    geomLocal = 0.d0
    R         = 0.d0
    h         = 0.d0
    triaArea  = 0.d0
    triaN1x   = 0.d0
    triaN1y   = 0.d0
    triaN2x   = 0.d0
    triaN2y   = 0.d0
    triaN3x   = 0.d0
    triaN3y   = 0.d0
!
! - Vectors for the triangle defined at mid_edge points
!
    X12 = (geom(4)+geom(13)-geom(1)-geom(10))/2.d0
    Y12 = (geom(5)+geom(14)-geom(2)-geom(11))/2.d0
    Z12 = (geom(6)+geom(15)-geom(3)-geom(12))/2.d0
    X13 = (geom(7)+geom(16)-geom(1)-geom(10))/2.d0
    Y13 = (geom(8)+geom(17)-geom(2)-geom(11))/2.d0
    Z13 = (geom(9)+geom(18)-geom(3)-geom(12))/2.d0
    !WRITE(6,*) 'Aretes: ', X12, Y12, Z12, X13, Y13, Z13
!
! - Normal to triangle
!
    Nx        = sqrt(X12*X12+Y12*Y12+Z12*Z12)
    triaNormX = X12/Nx
    triaNormY = Y12/Nx
    triaNormZ = Z12/Nx
    !WRITE(6,*) 'Nx: ', Nx
!
! - First tangent of triangle
!
    triaTau1X  = Y12*Z13-Z12*Y13
    triaTau1Y  = Z12*X13-X12*Z13
    triaTau1Z  = X12*Y13-Y12*X13
    Nz         = sqrt(triaTau1X*triaTau1X+triaTau1Y*triaTau1Y+triaTau1Z*triaTau1Z)
    !WRITE(6,*) 'Nz: ', Nz
    triaTau1X  = triaTau1X/Nz
    triaTau1Y  = triaTau1Y/Nz
    triaTau1Z  = triaTau1Z/Nz
!
! - Second tangent of triangle
!
    triaTau2X  = triaTau1Y*triaNormZ-triaTau1Z*triaNormY
    triaTau2Y  = triaTau1Z*triaNormX-triaTau1X*triaNormZ
    triaTau2Z  = triaTau1X*triaNormY-triaTau1Y*triaNormX
    Ny         = sqrt(triaTau2X*triaTau2X+triaTau2Y*triaTau2Y+triaTau2Z*triaTau2Z)
    !WRITE(6,*) 'Ny: ', Ny
    triaTau2X  = triaTau2X/Ny
    triaTau2Y  = triaTau2Y/Ny
    triaTau2Z  = triaTau2Z/Ny
!
! - Basis for mid_edge triangle
!
    R(1, 1:3) = (/triaNormX, triaNormY, triaNormZ/)
    R(2, 1:3) = (/triaTau2X, triaTau2Y, triaTau2Z/)
    R(3, 1:3) = (/triaTau1X, triaTau1Y, triaTau1Z/)
    !WRITE(6,*) 'R1: ',(/triaNormX, triaNormY, triaNormZ/)
    !WRITE(6,*) 'R2: ',(/triaTau2X, triaTau2Y, triaTau2Z/)
    !WRITE(6,*) 'R3: ',(/triaTau1X, triaTau1Y, triaTau1Z/)
!
! - Coordinates of element in the frame of mid_edge triangle
!
    geomLocal(1 ) = triaNormX*geom(1)+triaNormY*geom(2)+triaNormZ*geom(3)
    geomLocal(2 ) = triaTau2X*geom(1)+triaTau2Y*geom(2)+triaTau2Z*geom(3)
    geomLocal(3 ) = triaTau1X*geom(1)+triaTau1Y*geom(2)+triaTau1Z*geom(3)
    geomLocal(4 ) = triaNormX*geom(4)+triaNormY*geom(5)+triaNormZ*geom(6)
    geomLocal(5 ) = triaTau2X*geom(4)+triaTau2Y*geom(5)+triaTau2Z*geom(6)
    geomLocal(6 ) = triaTau1X*geom(4)+triaTau1Y*geom(5)+triaTau1Z*geom(6)
    geomLocal(7 ) = triaNormX*geom(7)+triaNormY*geom(8)+triaNormZ*geom(9)
    geomLocal(8 ) = triaTau2X*geom(7)+triaTau2Y*geom(8)+triaTau2Z*geom(9)
    geomLocal(9 ) = triaTau1X*geom(7)+triaTau1Y*geom(8)+triaTau1Z*geom(9)
    geomLocal(10) = triaNormX*geom(10)+triaNormY*geom(11)+triaNormZ*geom(12)
    geomLocal(11) = triaTau2X*geom(10)+triaTau2Y*geom(11)+triaTau2Z*geom(12)
    geomLocal(12) = triaTau1X*geom(10)+triaTau1Y*geom(11)+triaTau1Z*geom(12)
    geomLocal(13) = triaNormX*geom(13)+triaNormY*geom(14)+triaNormZ*geom(15)
    geomLocal(14) = triaTau2X*geom(13)+triaTau2Y*geom(14)+triaTau2Z*geom(15)
    geomLocal(15) = triaTau1X*geom(13)+triaTau1Y*geom(14)+triaTau1Z*geom(15)
    geomLocal(16) = triaNormX*geom(16)+triaNormY*geom(17)+triaNormZ*geom(18)
    geomLocal(17) = triaTau2X*geom(16)+triaTau2Y*geom(17)+triaTau2Z*geom(18)
    geomLocal(18) = triaTau1X*geom(16)+triaTau1Y*geom(17)+triaTau1Z*geom(18)
!
! - Height of element (thickness of shell)
!
    h = abs(geomLocal(12)+geomLocal(15)+geomLocal(18)&
           -geomLocal(3 )-geomLocal(6 )-geomLocal(9 ))/3.d0
    !WRITE(6,*) 'H: ',h
!
! - Coordinates of mid_edge triangle
!
    triaN1x = (geomLocal(1 )+geomLocal(10))/2.d0
    triaN1y = (geomLocal(2 )+geomLocal(11))/2.d0
    triaN2x = (geomLocal(4 )+geomLocal(13))/2.d0
    triaN2y = (geomLocal(5 )+geomLocal(14))/2.d0
    triaN3x = (geomLocal(7 )+geomLocal(16))/2.d0
    triaN3y = (geomLocal(8 )+geomLocal(17))/2.d0
    triaN3x = triaN3x - triaN1x
    triaN3y = triaN3y - triaN1y
    triaN2x = triaN2x - triaN1x
    triaN2y = 0.d0
    triaN1x = 0.d0
    triaN1y = 0.d0
    if (present(triaN1x_)) then
        triaN1x_ = triaN1x
        triaN1y_ = triaN1y
        triaN2x_ = triaN2x
        triaN2y_ = triaN2y
        triaN3x_ = triaN3x
        triaN3y_ = triaN3y
    endif
!
! - Area of mid_edge triangle
!
    triaArea = abs(triaN2x*triaN3y)/2.d0
    if (present(triaArea_)) then
        triaArea_ = triaArea
    endif
!
end subroutine
