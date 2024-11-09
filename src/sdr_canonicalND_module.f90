! Copyright (c) 2012-2013, 2017, 2019 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2012-2013 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2012 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2013 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
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
!> This module provides the routine to load canonical geometrical objects 
!! to growing array of geometries and their positions in growing array in
!! growing array of spatialObj
module sdr_canonicalND_module
  use env_module,                only: rk, labelLen
  use tem_canonicalND_module,    only: tem_canonicalND_type,                   &
    &                                  tem_load_canonicalND
  use tem_tools_module,          only: upper_to_lower
  use tem_transformation_module, only: tem_transformation_type
  use tem_point_module,          only: tem_point_type,                         &
    &                                  grw_pointArray_type, append
  use tem_line_module,           only: tem_line_type,                          &
    &                                  grw_lineArray_type, append
  use tem_triangle_module,       only: tem_triangle_type,                      &
    &                                  grw_triangleArray_type, append
  use tem_box_module,            only: tem_box_type, grw_boxArray_type,        &
    &                                  init, append
  use tem_plane_module,          only: tem_plane_type  

  use sdr_spatialObj_module,     only: grw_spatialObjArray_type,               &
    &                                  sdr_spatialObj_type,                    &
    &                                  point, line, triangle, box,             &
    &                                  append

  use flu_binding,               only: flu_State
  use aotus_module,              only: flu_State, aot_get_val,                 &
    &                                  aoterr_Fatal, aoterr_NonExistent,       &
    &                                  aoterr_WrongType
  use aot_table_module,          only: aot_table_open, aot_table_close,        &
    &                                  aot_table_length

  implicit none

  public :: sdr_load_canonicalND

contains

  ! *****************************************************************************
  !> This routine loads canonical geometrial objects like point, line, plane 
  !! and box and add them to the growing array of each primitive geometries
  !! and the position of this geometries are stored in the growing array
  !! of spatialObject
  subroutine sdr_load_canonicalND(pntArray, lineArray, triArray, boxArray, &
    &                             spaObjArray, attr_pos, transform,        &
    &                             conf, thandle)
    ! --------------------------------------------------------------------------!
    !> growing array of points
    type( grw_pointArray_type ), intent(inout) :: pntArray
    !> growing array of line
    type( grw_lineArray_type ), intent(inout) :: lineArray
    !> growing array of triangles
    type( grw_triangleArray_type ), intent(inout) :: triArray
    !> growing array of boxes
    type( grw_boxArray_type ), intent(inout) :: boxArray
    !> growing array of geometrical objects.
    type( grw_spatialObjArray_type), intent(inout) :: spaObjArray
    !> Position of the attribute to connect this object to.
    integer, intent(in) :: attr_pos
    !> transformation for spatial object
    type(tem_transformation_type), intent(in) :: transform 
    !> lua state
    type(flu_state) :: conf 
    integer, intent(in) :: thandle !< handle for canonical objects
    ! --------------------------------------------------------------------------!
    integer :: iCano
    type( tem_canonicalND_type ), allocatable :: canoND(:)
    ! --------------------------------------------------------------------------!

    ! Load canonical object definition
    call tem_load_canonicalND(canoND, transform, conf, thandle )

    do iCano = 1, size(canoND)
      select case(upper_to_lower(trim(canoND(iCano)%kind)))
      case('point')
        call append_CanoNDPointToSdrPoint(canopoint = canoND(iCano)%point, &
          &                               pntArray = pntArray,             &
          &                               spaObjArray = spaObjArray,       &
          &                               attr_pos = attr_pos              )

      case('line')
        call append_CanoNDLineToSdrLine(canoline    = canoND(iCano)%line, &
          &                             lineArray   = lineArray,          &
          &                             spaObjArray = spaObjArray,        &
          &                             attr_pos    = attr_pos            )

      case('plane')
        call append_CanoNDPlaneToTriangle(plane       = canoND(iCano)%plane, &
          &                               triArray    = triArray,            &
          &                               spaObjArray = spaObjArray,         &
          &                               attr_pos    = attr_pos             )

      case('box')
        if(canoND(iCano)%only_surface) then
          !hollow box
          call append_CanoNDBoxToTriangle(canobox     = canoND(iCano)%box, &
            &                             triArray    = triArray,          & 
            &                             spaObjArray = spaObjArray,       &
            &                             attr_pos    = attr_pos           )
        else  
          !solid box
          call append_CanoNDBoxToSdrBox(canobox     = canoND(iCano)%box, &
            &                           boxArray    = boxArray,          & 
            &                           spaObjArray = spaObjArray,       &
            &                           attr_pos    = attr_pos           )
        endif    

      end select
    end do
  end subroutine sdr_load_canonicalND
  ! *****************************************************************************


  ! *****************************************************************************
  !> This routine convert canonical point to seeder point array
  subroutine append_CanoNDPointToSdrPoint(canopoint, pntArray, spaObjArray, &
    &                                     attr_pos)
    ! --------------------------------------------------------------------------!
    !> canonical point object type
    type( tem_point_type ), intent(in) :: canopoint
    !> growing array of points
    type( grw_pointArray_type ), intent(inout) :: pntArray
    !> growing array of geometrical objects.
    type( grw_spatialObjArray_type), intent(inout) :: spaObjArray
    !> spatial object type
    integer, intent(in) :: attr_pos
    ! --------------------------------------------------------------------------!
    !local variables
    type(sdr_spatialObj_type) :: spatialObj 
    ! --------------------------------------------------------------------------!

    spatialObj%attribute_position = attr_pos
    spatialObj%geometry_primitive = point

    !append point to point list
    call append( pntArray, canopoint )

    spatialObj%primitive_position = pntArray%nVals

    !add spatialObj to spatialObj list
    call append( spaObjArray, spatialObj )

  end subroutine append_CanoNDPointToSdrPoint
  ! *****************************************************************************


  ! *****************************************************************************
  !> This routine convert canonical line to seeder line array
  subroutine append_CanoNDLineToSdrLine(canoline, lineArray, spaObjArray, &
    &                                   attr_pos)
    ! --------------------------------------------------------------------------!
    !> canonical geometry line object type
    type( tem_line_type ), intent(in) :: canoline
    !> growing array of line
    type( grw_lineArray_type ), intent(inout) :: lineArray
    !> growing array of geometrical objects.
    type( grw_spatialObjArray_type), intent(inout) :: spaObjArray
    !> spatial object type
    integer, intent(in) :: attr_pos
    ! --------------------------------------------------------------------------!
    !local variables
    type( sdr_spatialObj_type ) :: spatialObj 
    ! --------------------------------------------------------------------------!

    spatialObj%attribute_position = attr_pos
    spatialObj%geometry_primitive = line

    !append point to point list
    call append(lineArray, canoline)

    spatialObj%primitive_position = lineArray%nVals

    !add spatialObj to spatialObj list
    call append(spaObjArray, spatialObj)

  end subroutine append_CanoNDLineToSdrLine
  ! *****************************************************************************


  ! *****************************************************************************
  !> This routine convert box to planes and then planes to triangles
  !! and add triangle to growing array and add position of triangle into 
  !! spatialObj%primitive_position
  subroutine append_CanoNDBoxToTriangle(canoBox, triArray, spaObjArray, &
    &                                     attr_pos)
    ! --------------------------------------------------------------------------!
    !> canonical geometry box object type
    type( tem_box_type ), intent(in) :: canoBox
    !> growing array of triangles
    type( grw_triangleArray_type ), intent(inout) :: triArray
    !> growing array of geometrical objects.
    type( grw_spatialObjArray_type), intent(inout) :: spaObjArray
    !> spatial object type
    integer, intent(in) :: attr_pos
    ! --------------------------------------------------------------------------!
    !local variables
    ! --------------------------------------------------------------------------!

    !plane 1
    call append_CanoNDPlaneToTriangle(plane       = canoBox%plane(1), &
      &                               triArray    = triArray,         &
      &                               spaObjArray = spaObjArray,      &
      &                               attr_pos    = attr_pos          )

    !plane 2
    call append_CanoNDPlaneToTriangle(plane       = canoBox%plane(2), &
      &                               triArray    = triArray,         &
      &                               spaObjArray = spaObjArray,      &
      &                               attr_pos    = attr_pos          )

    !plane 3
    call append_CanoNDPlaneToTriangle(plane       = canoBox%plane(3), &
      &                               triArray    = triArray,         &
      &                               spaObjArray = spaObjArray,      &
      &                               attr_pos    = attr_pos          )

    !plane 4
    call append_CanoNDPlaneToTriangle(plane       = canoBox%plane(4), &
      &                               triArray    = triArray,         &
      &                               spaObjArray = spaObjArray,      &
      &                               attr_pos    = attr_pos          )

    !plane 5
    call append_CanoNDPlaneToTriangle(plane       = canoBox%plane(5), &
      &                               triArray    = triArray,         &
      &                               spaObjArray = spaObjArray,      &
      &                               attr_pos    = attr_pos          )

    !plane 6
    call append_CanoNDPlaneToTriangle(plane       = canoBox%plane(6), &
      &                               triArray    = triArray,         &
      &                               spaObjArray = spaObjArray,      &
      &                               attr_pos    = attr_pos          )

  end subroutine append_canoNDBoxToTriangle
  ! *****************************************************************************


  ! *****************************************************************************
  !> This routine convert plane to triangle and add triangle to 
  !! growing array and add position of triangle into 
  !! spatialObj%primitive_position
  subroutine append_CanoNDPlaneToTriangle(plane, triArray, spaObjArray, &
    &                                     attr_pos)
    ! --------------------------------------------------------------------------!
    !> canonical plane geometry object type
    type( tem_plane_type ), intent(in) :: plane
    !> growing array of triangles
    type( grw_triangleArray_type ), intent(inout) :: triArray
    !> growing array of geometrical objects.
    type( grw_spatialObjArray_type), intent(inout) :: spaObjArray
    !> spatial object type
    integer, intent(in) :: attr_pos
    ! --------------------------------------------------------------------------!
    !local variables
    type( sdr_spatialObj_type ) :: spatialObj 
    integer :: iTri
    ! --------------------------------------------------------------------------!
 
    spatialObj%attribute_position = attr_pos
    spatialObj%geometry_primitive = triangle

    !loop over triangles
    do iTri=1,2 ! two triangles per plane

      !append triangle into triangle list
      call append( triArray, plane%triangle(iTri) )

      spatialObj%primitive_position = triArray%nVals

      !add spatialObj to spatialObj list
      call append(spaObjArray, spatialObj)
    end do  

  end subroutine append_canoNDPlaneToTriangle
  ! *****************************************************************************


  ! *****************************************************************************
  !> This routine convert plane to triangle and add triangle to 
  !! growing array and add position of triangle into 
  !! spatialObj%primitive_position
  subroutine append_CanoNDBoxToSdrBox(canobox, boxArray, spaObjArray, &
    &                                 attr_pos)
    ! --------------------------------------------------------------------------!
    !> canonical geometry object type
    type( tem_box_type ), intent(in) :: canobox
    !> growing array of boxes
    type( grw_boxArray_type ), intent(inout) :: boxArray
    !> growing array of geometrical objects.
    type( grw_spatialObjArray_type), intent(inout) :: spaObjArray
    !> spatial object type
    integer, intent(in) :: attr_pos
    ! --------------------------------------------------------------------------!
    !local variables
    type( sdr_spatialObj_type ) :: spatialObj 
    ! --------------------------------------------------------------------------!
     
    spatialObj%attribute_position = attr_pos
    spatialObj%geometry_primitive = box

    call append(boxArray, canobox)

    spatialObj%primitive_position = boxArray%nVals

    ! Add spatialObj to spatialObj list.
    call append(spaObjArray, spatialObj)

  end subroutine append_canoNDBoxToSdrBox
  ! *****************************************************************************

end module sdr_canonicalND_module
