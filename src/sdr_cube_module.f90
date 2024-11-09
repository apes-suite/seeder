! Copyright (c) 2011-2014 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2011 Metin Cakircali <m.cakircali@grs-sim.de>
! Copyright (c) 2012-2013 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2012 Simon Zimny <s.zimny@grs-sim.de>
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
!> summary: This module contains cube definitions and routine to load
!! cube table from config file
! ******************************************************************************
module sdr_cube_module
  use env_module,         only: rk, labelLen
  use aotus_module,       only: flu_State, aot_get_val,                        &
    &                           aoterr_Fatal, aoterr_NonExistent,              &
    &                           aoterr_WrongType
  use aot_table_module,   only: aot_table_open, aot_table_close,               &
    &                           aot_table_length
  use tem_aux_module,     only: tem_abort
  use tem_logging_module, only: logunit

  implicit none

  !> An auxilary data type to describe a cube.
  !!
  !! Here: The origin is the corner from which on the cube is spanned 
  !! with the given length in each direction.
  type sdr_cube_type
    real(kind=rk) :: origin(3) !< origin of the cube
    real(kind=rk) :: center(3) !< center of the cube
    real(kind=rk) :: extent !< length of the cube
    real(kind=rk) :: halfwidth !< half length of the cube
  end type sdr_cube_type

contains

  ! ****************************************************************************
  !> This routine loads the boundCube table from config file
  ! ****************************************************************************
  subroutine sdr_load_cube(me, conf, key, pos, parent)
    ! --------------------------------------------------------------------------!
    type(sdr_cube_type), intent(out) :: me !< seeder cube type
    type(flu_state) :: conf !< lua state
    !> open cube table by given key
    character(len=*), optional, intent(in) :: key
    !> open cube table by position
    integer, optional, intent(in) :: pos
    !> if cube is to be load from pos, parent handle is required
    integer, optional, intent(in) :: parent
    ! --------------------------------------------------------------------------!
    integer :: thandle !< cube handle
    integer :: iError, vError(3), errFatal(3)
    real(kind=rk) :: deflen
    ! --------------------------------------------------------------------------!
    errFatal = aoterr_fatal

    me%origin = 0.0
    deflen = 1.0_rk
   
    if (present(key)) then
      ! open cube table from given key
      call aot_table_open(L=conf, thandle=thandle, key=trim(key))
      if (thandle == 0) then
        write(logunit(0),*) 'FATAL Error: cube definition not found with key:'
        write(logunit(0),*) (trim(key))
        call tem_abort()
      end if
    else if (present(parent) .and. present(pos)) then
      ! else open cube from given table pos
      call aot_table_open(L=conf, parent = parent, thandle=thandle, pos=pos)
      if (thandle == 0) then
        write(logunit(0),*) &
          &  'FATAL Error: cube definition not found with pos:', pos
        call tem_abort()
      end if
    else
      write(logunit(0),*) 'ERROR: Neither key nor pos is provided to load ' &
        &                 //' cube table'
      call tem_abort()
    end if

    call aot_get_val(L=conf, thandle=thandle, val=me%origin, ErrCode=vError, &
      &              key='origin', pos = 1, default=[0.0_rk, 0.0_rk, 0.0_rk] )
    if (any(btest(vError, errFatal))) then
      write(logunit(0),*) 'FATAL Error occured, while retrieving cube origin.'
      call tem_abort()
    end if

    call aot_get_val(L=conf, thandle=thandle, &
      &              val=me%extent, ErrCode=iError, &
      &              key='length', pos=2, default=deflen)
    if (btest(iError, aoterr_Fatal)) then
      write(logunit(0),*) 'FATAL Error occured, while retrieving cube length:'
      if (btest(iError, aoterr_NonExistent)) &
        &  write(logunit(0),*) 'Variable not existent!'
      if (btest(iError, aoterr_WrongType)) &
        &  write(logunit(0),*) 'Variable has wrong type!'
      call tem_abort()
    else
      if (btest(iError, aoterr_NonExistent)) &
        & write(logunit(1),*) 'Variable length not set in configuration, ' &
        &                     //'using default value!'
    end if

    me%halfwidth = 0.5_rk*me%extent
    me%center = me%origin + me%halfwidth

    write(logunit(1),*) ' Cube origin: ', me%origin
    write(logunit(1),*) ' Cube length: ', me%extent

    call aot_table_close(L=conf, thandle=thandle)

  end subroutine sdr_load_cube
  ! ****************************************************************************

end module sdr_cube_module
