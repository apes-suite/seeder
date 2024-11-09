! Copyright (c) 2011-2014, 2017 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2011-2013, 2019 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2011 Metin Cakircali <m.cakircali@grs-sim.de>
! Copyright (c) 2012 Sathish Krishnan P S <s.krishnan@grs-sim.de>
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
!> author: Kannan Masilamani
!! This module contain adapter routine to load sphere objects and add them to
!! list of spatial objects

module sdr_sphere_module
  use env_module,                only: rk
  use tem_aux_module,            only: tem_abort
  use tem_logging_module,        only: logunit
  use tem_transformation_module, only: tem_transformation_type  
  use tem_sphere_module,         only: tem_sphere_type, tem_load_sphere, &
    &                                  grw_sphereArray_type, append

  use aotus_module,              only: flu_State, aot_get_val, &
    &                                  aoterr_Fatal, aoterr_NonExistent,       &
    &                                  aoterr_WrongType
  use aot_table_module,          only: aot_table_open, aot_table_close,        &
    &                                  aot_table_length

  use sdr_spatialObj_module,     only: grw_spatialObjArray_type,               &
    &                                  sdr_spatialObj_type, append,            &
    &                                  sphere

  implicit none
  private

  public :: sdr_load_sphere

contains

  ! ****************************************************************************
  !> Load sphere information from config file.
  subroutine sdr_load_sphere(sphArray, spaObjArray, attr_pos, transform, &
    &                        conf, thandle )
    ! -------------------------------------------------------------------------!
    !inferface variables
    !> growing array of spheres
    type(grw_sphereArray_type), intent(inout) :: sphArray
    !> growing array of geometrical objects.
    type(grw_spatialObjArray_type), intent(inout) :: spaObjArray
    !> Position of the attribute to connect this object to.
    integer, intent(in) :: attr_pos
    !> transformation for spatial object
    type(tem_transformation_type), intent(in) :: transform 
    !> lua state
    type(flu_state) :: conf 
    integer, intent(in) :: thandle !< handle for canonical objects
    ! -------------------------------------------------------------------------!
    ! local varaibles
    integer :: iObj
    type(tem_sphere_type), allocatable :: loc_sphere(:)
    type( sdr_spatialObj_type ) :: spatialObj 
    ! -------------------------------------------------------------------------!

    ! load sphere objects
    call tem_load_sphere( loc_sphere, transform, conf, thandle )

    do iObj = 1, size(loc_sphere)
      !append sphere to sphere array
      spatialObj%attribute_position = attr_pos
      spatialObj%geometry_primitive = sphere
  
      call append( sphArray, loc_sphere(iObj) )
  
      spatialObj%primitive_position = sphArray%nVals
  
      !add spatialObj to spatialObj list
      call append(spaObjArray, spatialObj)
    end do

  end subroutine sdr_load_sphere
  ! ****************************************************************************

end module sdr_sphere_module

!> \page sphere Sphere
!! Spheres are defined by an origin and radius.
!! Sphere is considered to be solid as default i.e. all the cubes inside the
!! sphere are marked as intersected cubes. 
!! It is possible to created hollow spheres by setting only_surface = true,
!! it will mark only the cubes intersect with sphere surface as intersected
!! cubes
!!
!! Valid definition:
!! \li Single sphere
!! \verbatim
!! geometry={
!!   kind='sphere', 
!!     object={
!!       origin={0.0,0.0,0.0},
!!       radius=0.25,
!!       only_surface = true, -- If not defined default is set to false
!!     }
!! }
!! \endverbatim
!! 
!! \li Multiple sphere
!! \verbatim
!! geometry={
!!   kind='sphere', 
!!     object={
!!       {
!!       origin={0.0,0.0,0.0},
!!       radius=0.25
!!       },
!!       {
!!       origin={-2.0,0.0,0.0},
!!       radius=0.25
!!       }
!!     }
!! } 
!! \endverbatim
!! \n\n
!! The following seeder file is to generate mesh with hollow sphere (hollow => only_surface=true)
!!  inside:
!! \include testsuite/sphere/seeder.lua
!! \n\n
!! Mesh generated with hollow sphere by the seeder file:
!! \image html sphere.png
!! \image html sphere_withedges.png
!! \n\n
!! Cutview of mesh with hollow sphere:
!! \image html sphere_hollow.png
!! \n\n
!! Cutview of mesh with solid sphere (solid => only_surface=false):
!! \image html sphere_solid.png
!! \n\n
!! Example lua file is available at \link testsuite/sphere/seeder.lua
!! \example testsuite/sphere/seeder.lua

