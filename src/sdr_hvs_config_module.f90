! Copyright (c) 2015-2016 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2015 Samuel Ziegenberg
! Copyright (c) 2015 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2016 Neda Ebrahimi Pour <neda.epour@uni-siegen.de>
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
!> Harvesting configuration module for seeder.
!!
!! The configuration consists of the mesh to read and visualize, the output
!! definition providing the format to write, data extraction via trackings
!! and the configuration for the sampling of polynomials.
module sdr_hvs_config_module
  use env_module, only: PathLen

  use flu_binding, only: flu_State
  use aotus_module, only: aot_get_val, open_config_file, close_config, &
    &                     aoterr_Fatal

  use treelmesh_module, only: treelmesh_type, load_tem
  use tem_aux_module, only: tem_abort
  use tem_bc_prop_module, only: tem_bc_prop_type, init_tem_bc_prop
  use tem_color_prop_module, only: tem_color_prop_type, tem_color_prop_load
  use tem_general_module, only: tem_general_type, tem_load_general
  use tem_logging_module,   only: logunit, tem_logging_load_primary
  use tem_tracking_module, only: tem_tracking_type, tem_trackingControl_type, &
    &                            tem_load_tracking
  use tem_varsys_module, only: tem_varsys_type
  use tem_time_module,          only: tem_time_type, tem_time_reset
  use tem_restart_module,    only: tem_restart_type, tem_restart_openRead, &
    &                              tem_load_restart
  use tem_debug_module,     only: tem_debug_load_main
  use ply_sampled_tracking_module, only: ply_sampled_tracking_type, &
    &                                    ply_sampled_tracking_load

  use hvs_output_module, only: hvs_output_config_type, hvs_output_load

  use sdr_hvs_props_module, only: sdr_hvs_props_type, sdr_hvs_props_load

  implicit none

  !> This datatype describes the various settings to load from the configuration
  !! file.
  type sdr_hvs_config_type
    !> Location on disk to load the mesh data from.
    !!
    !! This prefix will be put before the various filenames of the individual
    !! mesh data files.
    character(len=pathLen) :: prefix

    !> Indication on wether to do subsampling.
    !!
    !! Subsampling is only necessary if there is any variable with more than a
    !! single degree of freedom, and the user requests sampling of the
    !! polynomials.
    logical :: do_subsampling

    !> Definition of trackings to extract only parts of the information.
    !!
    !! This includes the configuration for subsampling (which is only relevant
    !! if there is actually subresolved data in the mesh).
    type(ply_sampled_tracking_type) :: ply_sample_track

    !> Description of how the visualization output should be done.
    type(hvs_output_config_type) :: output

  end type sdr_hvs_config_type


contains


  !----------------------------------------------------------------------------!
  !> Read the configuration for the Seeder harvesting from a Lua script.
  subroutine sdr_hvs_config_load( me, mesh, property, varsys, general, &
    &                             restart, time                        )
  !----------------------------------------------------------------------------!
    !> Seeder harvesting configuration to load
    type(sdr_hvs_config_type), intent(out) :: me

    !> Treelm mesh description as obtained from the user configuration.
    type(treelmesh_type), intent(out) :: mesh

    !> Properties associated with the mesh.
    type(sdr_hvs_props_type), intent(out) :: property

    !> A variable system, to which the further variables should be appended.
    type(tem_varsys_type), intent(inout) :: varsys

    !> General treelm data to load
    type(tem_general_type), intent(inout) :: general

    !> 
    type(tem_restart_type), intent(inout) :: restart
 
    !>
    type(tem_time_type), intent(out)      :: time
    !----------------------------------------------------------------------!
    character(len=pathLen) :: filename
    type(flu_State) :: conf
    integer :: iError
    !----------------------------------------------------------------------!

    ! Load configuration data according to command line arguments.
    filename = ''

    ! Get filename from command argument
    call get_command_argument(1, filename)

    if ( trim(filename) == '' )  then
      ! Default to harvester.lua, if no filename is provided on the command
      ! line.
      filename = 'harvester.lua'
    end if

    ! Attempt to open the given file as a Lua script. Store a handle to that
    ! script in conf.
    call open_config_file(L = conf, filename = trim(filename))

    ! load and initialize logUnit
    call tem_logging_load_primary(conf = conf,             &
      &                           rank = general%proc%rank )

    ! load and initialize debug unit
    call tem_debug_load_main(conf = conf,                &
      &                      rank = general%proc%rank )

    call tem_load_general(me = general, conf = conf)

    call tem_load_restart( me       = restart,     &
      &                    conf     = conf,        &
      &                    tree     = mesh,        &
      &                    timing   = time,        &
      &                    globProc = general%proc )

    if (.not. restart%controller%readRestart) then
      call load_tem( me     = mesh,                   &
        &            conf   = conf,                   &
        &            myPart = general%proc%rank,      &
        &            nParts = general%proc%comm_size, &
        &            comm   = general%proc%comm       )
    end if

    call sdr_hvs_props_load( me     = property,    &
      &                      mesh   = mesh,        &
      &                      varsys = varsys,      &
      &                      proc   = general%proc )


    call ply_sampled_tracking_load( me   = me%ply_sample_track, &
      &                             conf = conf                 )

    ! Subresolved material data is the only information that provides details
    ! on a subelement resolution and might need sampling of polynomials.
    ! Activate the sampling only if there are actually subresolution information
    ! with non-constant data.
    me%do_subsampling = ( (me%ply_sample_track%sampling%max_nlevels > 0) &
      &                   .and. (property%subres%polydegree > 0) )

    if (me%do_subsampling) then
      write(logunit(1),*) 'Subresolved data will be sampled with at most ', &
        &                 me%ply_sample_track%sampling%max_nlevels, ' levels.'
    else
      me%ply_sample_track%sampling%max_nlevels = 0
      if (property%subres%polydegree > 0) then
        write(logunit(1),*) 'NO sampling of polynomial data will be done!'
      else
        write(logunit(3),*) 'NO polynomial data that would require subsampling.'
      end if
    end if

    ! If tracking table is not defined, load output_folder key
    ! to dump restart input or mesh to disk
    if (.not. me%ply_sample_track%tracking%control%active) then
      ! load output format and other config for output from output table
      call hvs_output_load( me       = me%output, &
        &                   conf     = conf,      &
        &                   isReduce = .false.    )

      ! Load output folder 
      call aot_get_val( L       = conf,            &
        &               key     = 'output_folder', &
        &               val     = me%prefix,       &
        &               default = './',            &
        &               ErrCode = iError           ) 
      write(logUnit(1),*) 'Ouput folder: '//trim(me%prefix)
 
    end if    

    ! Close the configuration script again.
    call close_config(conf)

    write(logunit(5),*) 'Done in sdr_hvs_config_load.'

  end subroutine sdr_hvs_config_load
  !----------------------------------------------------------------------------!
  !----------------------------------------------------------------------------!

end module sdr_hvs_config_module
