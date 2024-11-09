! Copyright (c) 2012-2013, 2019 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2012-2014 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2012 Sathish Krishnan P S <s.krishnan@grs-sim.de>
! Copyright (c) 2013 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
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
! ******************************************************************************
!> This module provides a adapter rouinte to load stl file and add triangles
!! to spatial objects 
!!
!! \author Kannan Masilamani
!!
module sdr_stl_module

  ! include treelm modules
  use env_module,                only: rk, PathLen
  use tem_aux_module,            only: tem_abort
  use tem_tools_module,          only: tem_horizontalSpacer
  use tem_logging_module,        only: logunit
  use tem_stlb_io_module,        only: tem_read_stlb, tem_size_stlb
  use tem_triangle_module,       only: tem_triangle_type,                      &
    &                                  grw_triangleArray_type, append
  use tem_transformation_module, only: tem_transformation_type  
  use tem_stl_module,            only: tem_load_stl, tem_stlData_type

  ! include seeder modules
  use stla_io,                   only: stla_read, stla_size, stla_check
  use sdr_spatialObj_module,     only: grw_spatialObjArray_type,               &
    &                                  sdr_spatialObj_type, triangle, append

  ! include aotus modules
  use flu_binding,               only: flu_State
  use aotus_module,              only: flu_State, aot_get_val,                 &
    &                                  aoterr_Fatal, aoterr_NonExistent,       &
    &                                  aoterr_WrongType
  use aot_table_module,          only: aot_table_open, aot_table_close,        &
    &                                  aot_table_length

  implicit none

  private

  public :: sdr_load_stl


contains

  ! *****************************************************************************
  !> This routine loads STL files from config and reads the triangles from the
  !! files into the dynamic array of triangles.
  subroutine sdr_load_stl(triArray, spaObjArray, attr_pos, transform,          &
    &                     conf, thandle)
    ! --------------------------------------------------------------------------!
    !> Dynamic array of triangles
    type(grw_triangleArray_type), intent(inout) :: triArray
    !> Growing array of geometrical objects.
    type(grw_spatialObjArray_type), intent(inout) :: spaObjArray
    !> Position of the attribute, this object is connected to.
    integer, intent(in) :: attr_pos
    !> transformation for spatial object
    type(tem_transformation_type), intent(in) :: transform 
    !> Lua state
    type(flu_state) :: conf 
    integer, intent(in) :: thandle !< handle for canonical objects
    ! --------------------------------------------------------------------------!
    integer :: iTri
    type(tem_stlData_type) :: stl_data
    type(sdr_spatialObj_type) :: spatialObj 
    type(tem_triangle_type) :: loc_tri  
    ! --------------------------------------------------------------------------!

    call tem_load_stl( stl_Data, transform, conf, thandle )

    ! Fix spatial object entries for all triangles from this STL.
    spatialObj%attribute_position = attr_pos
    spatialObj%geometry_primitive = triangle

    do iTri=1,stl_data%nTris
      loc_tri%nodes(:,1) = stl_data%nodes( :, stl_data%tri_node(1,iTri) )
      loc_tri%nodes(:,2) = stl_data%nodes( :, stl_data%tri_node(2,iTri) )
      loc_tri%nodes(:,3) = stl_data%nodes( :, stl_data%tri_node(3,iTri) )

      ! Append triangle into list of all triangles.
      call append(triArray, loc_tri)

      ! Store the position of this triangle for identification in the list of
      ! all spatial objectes.
      spatialObj%primitive_position = triArray%nVals

      ! Append spatialObj to list of spatialObjs.
      call append(spaObjArray, spatialObj)
    end do

    ! Deallocate nodes and triangle array in stl_data.
    deallocate(stl_data%nodes)
    deallocate(stl_data%tri_node)

  end subroutine sdr_load_stl
  ! *****************************************************************************

end module sdr_stl_module

!> \page stl STL
!! stl files can be used as geometry kind. At first, seeder load the 
!! triangles from stl files to the temporary stl_data type and
!! then each triangle in the stl_data is converted to sdr_triangle_type.\n
!! stl geometry requires filename and stl file format. If file format
!! is not provided, default is set to binary.
!! Valid definition:
!! \li Single stl
!! \verbatim
!! geometry={
!!   kind='stl', 
!!     object={
!!       filename='cube.stl',
!!       format = 'ascii' -- if not provided, default is binary
!!     }
!! }
!! \endverbatim
!! 
!! \li Multiple stls
!! \verbatim
!! geometry={
!!   kind='stl', 
!!     object={
!!       {
!!       filename = 'cube.stl'
!!       },
!!       {
!!       filename = 'cylinder.stl'
!!       }
!!     }
!! }
!! \endverbatim 
!! \n\n
!! Seeder file to create mesh with single 'stl' geometry:
!! \include testsuite/stl/seeder.lua
!! \n\n
!! Mesh with 'stl' geometry created by seeder file:
!! \image html stl_cylinder.png
!! \endverbatim
!! Example lua file is available at \link testsuite/stl/seeder.lua
!! \example testsuite/stl/seeder.lua
 
