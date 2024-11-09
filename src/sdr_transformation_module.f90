! Copyright (c) 2012-2013 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2013-2014, 2017 Harald Klimach <harald.klimach@uni-siegen.de>
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
! *****************************************************************************!
!> This module provide datatype and routine for transformation of
!! geometrical objects i.e translation and deformation
module sdr_transformation_module
  use env_module,             only: rk
  use tem_aux_module,         only: tem_abort
  use tem_tools_module,       only: upper_to_lower
  use tem_logging_module,     only: logunit
  use tem_canonicalND_module, only: tem_canonicalND_type

  use flu_binding,            only: flu_State
  use aotus_module,           only: aot_get_val, flu_State,           &
    &                               aoterr_Fatal, aoterr_NonExistent, &
    &                               aoterr_WrongType
  use aot_table_module,       only: aot_table_open, aot_table_close, &
    &                               aot_table_length, aot_table_top, &
    &                               aot_table_first

  implicit none

  private 

  public :: sdr_transformation_type
  public :: sdr_translation_type
  public :: sdr_deformation_type
  public :: sdr_transformCanoND
  public :: sdr_load_transformation

  !> Data type defines geometry translation
  type sdr_translation_type
    !>Is translation defined
    logical :: active
    !> vector defining translation in x,y,z direction
    real(kind=rk) :: vec(3)
  end type sdr_translation_type

  !> Data type defines geometry scale and rotation
  type sdr_deformation_type
    !> Is deformation defined
    logical :: active
    !> matrix defining the deformation
    real(kind=rk) :: matrix(3,3)
  end type sdr_deformation_type

  !> Data type defines geometry transformation 
  type sdr_transformation_type
    !> is transformation active
    logical :: active
    !> translation of geometry
    type(sdr_translation_type) :: translate
    !> deformation of geometry
    type(sdr_deformation_type) :: deform
  end type sdr_transformation_type

  !> This routine apply transformations to canonical objects
  interface sdr_transformcanoND
    module procedure transformcanoND
    module procedure transformcanoND_single
  end interface sdr_transformcanoND


contains


  ! ************************************************************************** !
  !> This routine loads the transformation table for each spatial object table
  !! in config file
  !!
  !! If single spatial object contains multiple geometry then the transformation
  !! is applied to all geometries defined in that spatial object
  subroutine sdr_load_transformation( transform, conf, thandle )
    ! ---------------------------------------------------------------------- !
    !inferface variables
    !> transformation for spatial object
    type(sdr_transformation_type), intent(out) :: transform
    !> lua state
    type(flu_state) :: conf
    !> spatial object parent handle
    integer, intent(in) :: thandle
    ! ---------------------------------------------------------------------- !
    integer :: transform_handle
    ! ---------------------------------------------------------------------- !

    !set default to false
    transform%active = .false.

    call aot_table_open(L = conf, parent = thandle, &
      &                 thandle = transform_handle, &
      &                 key = 'transformation'      )


    if (transform_handle > 0) then
      write(logunit(1),*) 'Loading transformation '
      transform%active = .true.
      !load translation table
      call sdr_load_translation( translate = transform%translate, &
        &                        conf      = conf,                &
        &                        thandle   = transform_handle     )

      !load deformation table
      call sdr_load_deformation( deform  = transform%deform, &
        &                        conf    = conf,             &
        &                        thandle = transform_handle  )
     endif

    call aot_table_close(L=conf, thandle=transform_Handle)

  end subroutine sdr_load_transformation
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> This routine loads the translation table from transformation table
  subroutine sdr_load_translation( translate, conf, thandle )
    ! ---------------------------------------------------------------------- !
    !inferface variables
    !> translate for spatial object
    type(sdr_translation_type), intent(out) :: translate
    !> lua state
    type(flu_state) :: conf
    !> spatial object parent handle
    integer, intent(in) :: thandle
    ! ---------------------------------------------------------------------- !
    integer :: translate_handle
    integer :: vError(3), errFatal(3)
    ! ---------------------------------------------------------------------- !
    errFatal = aoterr_fatal

    translate%active = .false.
    translate%vec = 0.0_rk
    call aot_table_open(L = conf, parent = thandle, &
      &                 thandle = translate_handle, &
      &                 key = 'translation')

    if(translate_handle > 0) then
      translate%active = .true.

      call aot_get_val(L=conf, thandle = thandle, key = 'translation', &
        &              val=translate%vec, ErrCode = vError             )

      if (any(btest(vError, errFatal))) then
        write(logunit(0),*) 'Error in configuration: translate table in'
        write(logunit(0),*) '    transformation table'
        call tem_abort()
      end if
    endif

    call aot_table_close(L=conf, thandle=translate_Handle)

    if (translate%active) then
      write(logunit(1),*) ' Translation = ', translate%vec
    endif

  end subroutine sdr_load_translation
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> This routine loads the deformation table from transformation table
  subroutine sdr_load_deformation( deform, conf, thandle )
    ! ---------------------------------------------------------------------- !
    !inferface variables
    !> deform for spatial object
    type(sdr_deformation_type), intent(out) :: deform
    !> lua state
    type(flu_state) :: conf
    !> spatial object parent handle
    integer, intent(in) :: thandle
    ! ---------------------------------------------------------------------- !
    integer :: deform_handle, deform_subhandle
    integer :: iError, iPos
    integer :: vError(3), errFatal(3)
    real(kind=rk) :: const, vec(3)
    ! ---------------------------------------------------------------------- !
    errFatal = aoterr_fatal

    deform%active = .false.
    deform%matrix = 0.0_rk


    !First check if deformation is defined as a constant
    call aot_get_val(L=conf, thandle = thandle, key='deformation', &
      &              val=const, ErrCode=iError                     )
    if (btest(iError, aoterr_fatal)) then
      call aot_table_open(L = conf, parent = thandle, &
        &                 thandle = deform_handle,    &
        &                 key = 'deformation'         )

      if (deform_handle > 0) then
        deform%active = .true.
        !deformation is defined as a table
        call aot_table_open(L = conf, parent = deform_handle, &
          &                 thandle = deform_subhandle, &
          &                 pos = 1)
        if(deform_subhandle > 0) then
          call aot_table_close(L=conf, thandle=deform_subHandle)
          !sub table exist load it as matrix
          do iPos=1,3
            call aot_get_val(L=conf, thandle=deform_handle, &
              &              pos=iPos, val=vec, ErrCode=vError) 

            if (any(btest(vError, errFatal))) then
              write(logunit(0),*) 'Error in configuration: '
              write(logunit(0),*) '     Deformation table at pos', iPos
              call tem_abort()
            end if
            deform%matrix(iPos,:) = vec 
          end do
        else
          call aot_table_close(L=conf, thandle=deform_Handle)
          !if single table then it is vec
          call aot_get_val(L=conf, thandle=thandle, &
            &              key='deformation', val=vec, ErrCode=vError) 
          if (any(btest(vError, errFatal))) then
            write(logunit(0),*) 'Error in configuration: '
            write(logunit(0),*) '   Loading deformation table as vector'
            call tem_abort()
          end if
          deform%matrix(1,1) = vec(1)
          deform%matrix(2,2) = vec(2)
          deform%matrix(3,3) = vec(3)
        endif
      else !deformation table not defined
        !close the table
        call aot_table_close(L=conf, thandle=deform_Handle)
      endif
    else
      !it is a constant scaling factor
      deform%active = .true.
      deform%matrix(1,1) = const
      deform%matrix(2,2) = const
      deform%matrix(3,3) = const
    endif


    if (deform%active) then
      write(logunit(1),"(A,3E12.5)") ' Deformation = ', deform%matrix(1,:)
      write(logunit(1),"(A,3E12.5)") '               ', deform%matrix(2,:)
      write(logunit(1),"(A,3E12.5)") '               ', deform%matrix(3,:)
    endif

  end subroutine sdr_load_deformation
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> This routine apply transformation to canonical objects.
  subroutine transformCanoND(canoND, transform)
    ! ---------------------------------------------------------------------- !
    !> canonical geometry object type
    type( tem_canonicalND_type ), intent(inout) :: canoND(:)
    !> transformation for spatial object
    type(sdr_transformation_type), intent(in) :: transform
    ! ---------------------------------------------------------------------- !
    integer :: iCano
    ! ---------------------------------------------------------------------- !

    do iCano=1,size(canoND)
      call transformCanoND_single(canoND = canoND(iCano), &
        &                         transform = transform   )
    end do

  end subroutine transformCanoND
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> This routine apply transformation to canonical objects.
  subroutine transformCanoND_single(canoND, transform)
    ! ---------------------------------------------------------------------- !
    !> canonical geometry object type
    type( tem_canonicalND_type ), intent(inout) :: canoND
    !> transformation for spatial object
    type(sdr_transformation_type), intent(in) :: transform 
    ! ---------------------------------------------------------------------- !
    integer :: iVec
    ! ---------------------------------------------------------------------- !
    if(transform%active) then
      if(transform%deform%active) then
        canoND%origin = matmul( transform%deform%matrix, &
          &                     canoND%origin            )
        do iVec=1,canoND%nDim
          canoND%vec(:,iVec) = matmul( transform%deform%matrix, &
            &                          canoND%vec(:,iVec)       )
        enddo
      endif
      if(transform%translate%active) then
        canoND%origin = transform%translate%vec + canoND%origin
      endif
    endif

  end subroutine transformCanoND_single
  ! ************************************************************************** !


end module sdr_transformation_module

!> \page transformation Transformation
!! Transformation is used to scale, translate, rotate and reflect
!! the geometrical objects. Trasformation table is defined in the spatial
!! object table in the lua config file. If the geometry in the spatial
!! object contains multiple geometries then the transformation defined in that
!! spatial object is applied to all the geometries.\n
!! If both translation and deformation are defined for the geometry object
!! then the deformation is applied first and then the deformed geometry
!! is then translated.\n
!!
!! \li Translation
!! Translation is a table with three entries defining x,y,z coordinate values
!! to translate the geometrical object.
!! Gometry is translated just by adding the position of the geometry with 
!! given translation vector.\n
!! Example:
!! \verbatim
!! spatial_object={
!!   ...<attribute>...
!!   ...<geometry>...
!!   transformation={
!!     translation={0.0,2.0,0.0} -- translating the object along y-axis by 2.0
!!   }
!! }
!! \endverbatim
!! \n
!! \li Deformation
!! Deformation table can be used to scale, rotate and reflect the geometry.
!! Deformation cane be defined as const, vector and matrix. In the code,
!! it is converted to matrix with 3x3. Matrix is multiplied with a geometry
!! vector to scale, rotate or reflect depends on the given matrix.
!!
!! * Scaling
!! If deformation is const, then geometry is scaled in all three directions
!! with const and it is vector with three entries, then geometry is scaled in
!! x,y,z directions with different factor.\n
!! Example:
!! ** Constant scaling in three direction
!!\verbatim
!! spatial_object={
!!   ...<attribute>...
!!   ...<geometry>...
!!   transformation={
!!     deformation = 2.0, -- scaled in all direction by 2.0
!!   }
!!}
!!\endverbatim
!! ** Different scaling in three direction
!!\verbatim
!! spatial_object={
!!   ...<attribute>...
!!   ...<geometry>...
!!   transformation={
!!     deformation = {0.5,2.0,1.5}
!!   }
!!}
!!\endverbatim
!! * Reflection
!! Below example reflect the geometry object in y-axis
!!\verbatim
!! spatial_object={
!!   ...<attribute>...
!!   ...<geometry>...
!!   transformation={
!!     deformation = {1.0,-1.0,1.0}
!!   }
!!}
!!\endverbatim
!! *Rotation
!! Rotation is defined by the deformation table with 3x3 entries.
!! Below example rotate the geometry object in z-axis in anti-clockwise
!! direction by 45 degrees.
!!\verbatim
!! spatial_object={
!!   ...<attribute>...
!!   ...<geometry>...
!!   transformation={
!!     deformation = {
!!                    { 0.5*math.cos(45*math.pi/180), 
!!                      -0.5*math.sin(45*math.pi/180), 
!!                      0.0 },
!!                    { 0.5*math.sin(45*math.pi/180), 
!!                      0.5*math.cos(45*math.pi/180), 
!!                      0.0 },
!!                    { 0.0, 0.0, 0.5 }
!!     }
!!   }
!!}
!!\endverbatim
!! More information on rotatation matrix can be found in
!! http://en.wikipedia.org/wiki/Rotation_(mathematics)
!! \n
!! It is also possible to combine scaling, reflection and rotation in the 
!! deformation matrix.
!! Example lua file is available at \link testsuite/transform/seeder.lua
!! \example testsuite/transform/seeder.lua
