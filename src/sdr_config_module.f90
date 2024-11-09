! Copyright (c) 2012-2013 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2012-2014, 2016 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2012 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2012, 2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2015 Samuel Ziegenberg
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2016 Peter Vitt <peter.vitt2@uni-siegen.de>
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
!> This module provides the configuration of the mesh generation.
module sdr_config_module
  use flu_binding,          only: flu_State, flu_next
  use aotus_module,         only: aot_get_val, open_config_file, close_config, &
    &                             flu_State,                                   &
    &                             aoterr_Fatal, aoterr_NonExistent,            &
    &                             aoterr_WrongType
  use aot_table_module,     only: aot_table_open, aot_table_close, &
    &                             aot_table_length, aot_table_top,             &
    &                             aot_table_first

  use env_module,           only: PathLen, labelLen
  use tem_general_module,   only: tem_general_type, tem_load_general
  use tem_logging_module,   only: logunit, tem_logging_load_primary
  use tem_aux_module,       only: tem_abort
  use tem_tools_module,     only: upper_to_lower
  use tem_comm_module,      only: tem_commpattern_type
  use tem_dyn_array_module, only: dyn_labelArray_type, init, append, &
    &                             PositionOfVal
  use tem_debug_module,     only: tem_debug_load_main

  use sdr_geometry_module,  only: sdr_geometry_type, sdr_load_geometry
  use sdr_attribute_module, only: sdr_seed_object
  use sdr_subresolution_module, only: sdr_subresolution_type, &
    &                                 sdr_subresolution_load, &
    &                                 sdr_subresolution_encolor

  use sdr_timer_module,     only: timer_handle_loadconfig
  use tem_timer_module,     only: tem_startTimer,tem_stopTimer
 
  implicit none

  private

  public :: sdr_load_config
  public :: sdr_confHead_type

  !> This type contains basic information loaded from the config file.
  type sdr_confHead_type
    character(len=PathLen) :: comment    !< Comment line in the file
    character(len=PathLen) :: folder     !< folder where to save the mesh
    integer :: minlevel !< minumum refinement level of the fluid tree
    type(tem_general_type) :: general

    !> Information on resolution of boundaries within elements using
    !! polynomials.
    type(sdr_subresolution_type) :: subresolution

    !> List of colors to invert.
    !!
    !! Colors that are listed here, will get their definition inverted, such
    !! that filled areas are turned void and void areas will get filled if
    !! they are flooded by the none color. Inverting the none color is not
    !! possible.
    type(dyn_labelArray_type) :: inverted_colors
  end type sdr_confHead_type


contains


  ! *****************************************************************************
  !> Load the configuration from the Lua script provided on the command line
  !! or from seeder.lua by default, if no file name is given as program
  !! argument.
  !!
  !! The configuration needs to describe some general properties, like the
  !! directory where the resulting mesh is to be stored and the geometries,
  !! that are to be used to define the mesh.
  subroutine sdr_load_config( me, geometry )
    ! --------------------------------------------------------------------------!
    !> contains basic information from config file
    !  @todo: me was filled before so intent has to be inout
    type( sdr_confHead_type ), intent(inout) :: me
    !> contains all geometry object defined in the config file
    type( sdr_geometry_type ), intent(out) :: geometry
    ! --------------------------------------------------------------------------!
    ! local variables
    character(len=1024) :: filename
    type(flu_State) :: conf
    integer :: iError
    ! --------------------------------------------------------------------------!

    call tem_startTimer(timerHandle = timer_handle_loadconfig) 
 
    filename = ''
    ! Get filename from command argument
    call get_command_argument(1,filename)
    if ( trim(filename) == '')  then
      ! Default to seeder.lua, if no filename is provided on the command line.
      filename = 'seeder.lua'
    end if

    ! Attempt to open the given file as a Lua script. Store a handle to that
    ! script in conf.
    call open_config_file(L = conf, filename = trim(filename))

    ! load and initialize logUnit
    call tem_logging_load_primary(conf = conf,                &
      &                           rank = me%general%proc%rank )

    ! load and initialize debug unit
    call tem_debug_load_main(conf = conf,                &
      &                      rank = me%general%proc%rank )

    call tem_load_general(me = me%general, conf = conf)

    call aot_get_val(L=conf, key='comment', &
      &              val=me%comment, ErrCode = iError, default='')
    write(logunit(3),*) 'Comment in the file: '// trim(me%comment)

    ! Get the directory, where the resulting mesh is to be stored in.
    call aot_get_val(L=conf, key='folder', &
      &              val=me%folder, ErrCode = iError, default='mesh_')
    if (btest(iError, aoterr_WrongType)) then
      write(logunit(0),*) "ERROR: 'folder' to load mesh has wrong type!"
      call tem_abort()
    endif

    ! Load the minimum refinement level for the mesh
    call aot_get_val(L=conf, key='minlevel', &
      &              val=me%minlevel, ErrCode = iError, default=0)
    if (btest(iError, aoterr_WrongType)) then
      write(logunit(0),*) "ERROR: minLevel has wrong type!"
      call tem_abort()
    endif

    write(logunit(1),*) 'Mesh folder: '//trim(me%folder)
    write(logunit(2),*) 'minlevel: ', me%minlevel

    call sdr_subresolution_load( subres = me%subresolution, &
      &                          conf   = conf              )

    call sdr_load_inverted(me%inverted_colors, conf)

    ! Load geometry object from config file.
    call sdr_load_geometry( me            = geometry,                        &
      &                     subres_colors = me%subresolution%default_colors, &
      &                     invert_colors = me%inverted_colors,              &
      &                     conf          = conf                             )

    ! Close the configuration script again.
    call close_config(conf)

    ! Now go through all colors, and ensure, that there is a definition of
    ! values for each. If there is no definition given by the user, we default
    ! to a filling value of 1 and a void value of 0.
    write(logunit(5),*) 'Getting encoloring values for all colors'
    call sdr_subresolution_encolor( me     = me%subresolution,          &
      &                             colors = geometry                   &
      &                                      %attribute                 &
      &                                      %uni_name(sdr_seed_object) )

    write(logunit(5),*) 'Done in sdr_load_config'

    call tem_stopTimer( timerHandle = timer_handle_loadconfig )

  end subroutine sdr_load_config


  subroutine sdr_load_inverted(invlist, conf)
    !> List of labels for inverted colors.
    type(dyn_labelArray_type), intent(out) :: invlist
    !> Lua script handle to read the inversion list from.
    type(flu_State) :: conf

    character(labelLen) :: cLabel
    integer :: thandle
    integer :: nInverted
    integer :: iInv
    integer :: pos
    integer :: iError
    logical :: wasAdded

    call aot_table_open(L = conf, thandle = thandle, &
      &                 key = 'inverted_colors')

    if (thandle > 0) then

      nInverted = aot_table_length(L = conf, thandle = thandle)
      call init(invlist, length = nInverted)

      do iInv=1,nInverted
        call aot_get_val( L       = conf,    &
          &               thandle = thandle, &
          &               pos     = iInv,    &
          &               val     = clabel,  &
          &               ErrCode = iError   )
        if (btest(iError, aoterr_Fatal)) then
          write(logunit(0),*) &
            &  'FATAL Error occured while retrieving color name in the' &
            &  //' inversion table!'
          call tem_abort()
        end if
        call append( invlist, upper_to_lower(trim(clabel)), &
          &          pos = pos, wasAdded = wasAdded         )
      end do

    end if

    call aot_table_close(L = conf, thandle = thandle)

  end subroutine sdr_load_inverted

end module sdr_config_module
