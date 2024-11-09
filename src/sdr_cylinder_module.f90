! Copyright (c) 2012-2014, 2019 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2012-2014, 2022 Harald Klimach <harald.klimach@dlr.de>
! Copyright (c) 2012 Sathish Krishnan P S <s.krishnan@grs-sim.de>
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
! ******************************************************************************
!> author: Kannan Masilamani
!! This module contain adapter routine to load cylinder objects and add them to
!! list of spatial objects
module sdr_cylinder_module
  use env_module,                only: rk
  use tem_param_module,          only: PI
  use tem_aux_module,            only: tem_abort
  use tem_logging_module,        only: logunit
  use tem_transformation_module, only: tem_transformation_type  
  use tem_cylinder_module,       only: tem_cylinder_type, tem_load_cylinder, &
    &                                  grw_cylinderArray_type, append 

  use aotus_module,              only: flu_State, aot_get_val,                 &
    &                                  aoterr_Fatal, aoterr_NonExistent,       &
    &                                  aoterr_WrongType
  use aot_table_module,          only: aot_table_open, aot_table_close,        &
    &                                  aot_table_length

  use sdr_spatialObj_module,     only: grw_spatialObjArray_type,               &
    &                                  sdr_spatialObj_type, append,            &
    &                                  cylinder

  implicit none
  private

  public :: sdr_load_cylinder
  public :: append_cylinder2SpaObj

contains
  ! ****************************************************************************
  !> \brief Loading cylinder information from config file \n
  subroutine sdr_load_cylinder(cylArray, spaObjArray, attr_pos, transform,     &
    &                          conf, thandle )
    ! --------------------------------------------------------------------------!
    !inferface variables
    !> growing array of cylinders
    type(grw_cylinderArray_type), intent(inout) :: cylArray
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
    ! local varaibles
    integer :: iObj
    type(tem_cylinder_type), allocatable :: loc_cylinder(:)
    ! --------------------------------------------------------------------------!

    ! load cylinder
    call tem_load_cylinder(loc_cylinder, transform, conf, thandle)

    !append cylinder
    do iObj = 1, size(loc_cylinder)
      call append_cylinder2SpaObj( cylArray, spaObjArray, attr_pos, &
        &                          loc_cylinder(iObj)               )
    end do

  end subroutine sdr_load_cylinder
  ! *****************************************************************************


  ! *****************************************************************************
  !> This routine single cylinder from object table
  subroutine append_cylinder2SpaObj(cylArray, spaObjArray, attr_pos, &
    &                               loc_cylinder)
    ! --------------------------------------------------------------------------!
    !inferface variables
    !> growing array of cylinders
    type(grw_cylinderArray_type), intent(inout) :: cylArray
    !> growing array of geometrical objects.
    type( grw_spatialObjArray_type), intent(inout) :: spaObjArray
    !> Position of the attribute to connect this object to.
    integer, intent(in) :: attr_pos
    type(tem_cylinder_type), intent(in) :: loc_cylinder
    ! --------------------------------------------------------------------------!
    type( sdr_spatialObj_type ) :: spatialObj 
    ! --------------------------------------------------------------------------!
    !append cylinder to cylinder array
    spatialObj%attribute_position = attr_pos
    spatialObj%geometry_primitive = cylinder

    call append( cylArray, loc_cylinder )

    spatialObj%primitive_position = cylArray%nVals

    !add spatialObj to spatialObj list
    call append(spaObjArray, spatialObj)

  end subroutine append_cylinder2SpaObj

end module sdr_cylinder_module

!> \page cylinder Cylinders
!! Cylinders are defined by an origin, vector defining the length and the
!! axis and the radius.
!! Cylinder is considered to be solid as default i.e. all the cubes inside the
!! cylinder are marked as intersected cubes. 
!! It is possible to created hollow cylinders by setting only_surface = true,
!! it will mark only the cubes intersect with cylinder surface as intersected
!! cubes
!!
!! Valid definition:
!! \li Single cylinder
!! \verbatim
!! geometry={
!!   kind='cylinder', 
!!     object={
!!       origin={0.0,0.0,0.0},
!!       vec={1.0,0.0,0.0},
!!       radius=0.25,
!!       only_surface = true, -- If not defined default is set to false
!!     }
!! }
!! \endverbatim
!! 
!! \li Multiple cylinder
!! \verbatim
!! geometry={
!!   kind='cylinder', 
!!     object={
!!       {
!!       origin={0.0,0.0,0.0},
!!       vec={1.0,0.0,0.0},
!!       radius=0.25
!!       },
!!       {
!!       origin={0.0,0.0,0.0},
!!       vec={1.0,1.0,0.0},
!!       radius=0.25
!!       }
!!     }
!! } 
!! \endverbatim
!! \n\n 
!! Seeder file to generate mesh with single cylinder (only_surface=true) is below:
!! include testsuite/plane/seeder.lua
!! \n\n
!! Mesh with hollow cylinder (Hollow => only_surface = true)
!! \image html cylinder.png
!! \n\n
!! \image html cylinder_withedges.png
!! \n\n
!! Cutview of mesh with hollow cylinder
!! \image html cylinder_hollow.png
!! \n\n
!! As said earlier, cylinder can be created as solid one using 'only_surface=false'.
!! Cutview of Mesh generated with 'only_surface=false':
!! \image html cylinder_solid.png
!! \n\n
!! Example lua file is available at \link testsuite/cylinder/seeder.lua
!! \example testsuite/cylinder/seeder.lua

