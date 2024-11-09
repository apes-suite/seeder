! Copyright (c) 2012-2013, 2019 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2012-2014, 2017 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2012 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2012 Sathish Krishnan P S <s.krishnan@grs-sim.de>
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright notice, this
! list of conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright notice,
! this list of conditions and the following disclaimer in the documentation
! and/or other materials provided with the distribution.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
! DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
! FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
! SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
! OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!> \brief This module is an adapter to tem_triangle_module, Contains 
!! a routine to load triangle and added it to list of spatial objects
module sdr_triangle_module
  use tem_logging_module,        only: logunit
  use tem_transformation_module, only: tem_transformation_type  
  use tem_triangle_module,       only: tem_triangle_type,      &
    &                                  grw_triangleArray_type, &
    &                                  tem_load_triangle, append

  use aotus_module,              only: flu_State

  use sdr_spatialObj_module,     only: grw_spatialObjArray_type,               &
    &                                  sdr_spatialObj_type, append,            &
    &                                  triangle

  implicit none

  private

  public :: sdr_load_triangle

contains

  ! ****************************************************************************
  !> Load triangle information from config file.
  subroutine sdr_load_triangle(triArray, spaObjArray, attr_pos, transform, &
    &                          conf, thandle )
    !--------------------------------------------------------------------------!
    !inferface variables
    !> growing array of triangles
    type(grw_triangleArray_type), intent(inout) :: triArray
    !> growing array of geometrical objects.
    type( grw_spatialObjArray_type), intent(inout) :: spaObjArray
    !> Position of the attribute to connect this object to.
    integer, intent(in) :: attr_pos
    !> transformation for spatial object
    type(tem_transformation_type), intent(in) :: transform 
    !> lua state
    type(flu_state) :: conf 
    integer, intent(in) :: thandle !< handle for canonical objects
    !--------------------------------------------------------------------------!
    ! local varaibles
    type(tem_triangle_type), allocatable :: loc_triangle(:)
    type(sdr_spatialObj_type) :: spatialObj 
    integer :: iTri
    !--------------------------------------------------------------------------!

    call tem_load_triangle( me        = loc_triangle, &
      &                     transform = transform,    &
      &                     conf      = conf,         &
      &                     thandle   = thandle       )

    !append triangle to triangle array
    spatialObj%attribute_position = attr_pos
    spatialObj%geometry_primitive = triangle

    do iTri = 1, size(loc_triangle)
      call append( triArray, loc_triangle(iTri) )

      spatialObj%primitive_position = triArray%nVals

      !add spatialObj to spatialObj list
      call append(spaObjArray, spatialObj)
    end do

  end subroutine sdr_load_triangle
  ! ****************************************************************************

end module sdr_triangle_module
