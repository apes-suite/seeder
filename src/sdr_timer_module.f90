! Copyright (c) 2015 Samuel Ziegenberg
! Copyright (c) 2015 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2015-2016, 2022 Harald Klimach <harald.klimach@dlr.de>
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
module sdr_timer_module
  use tem_timer_module,        only: tem_addTimer, tem_getMaxTimerVal, &
    &                                tem_getTimerName
  use soi_revision_module,     only: soi_solver_revision
  use tem_general_module,      only: tem_general_type
  use env_module,              only: rk, pathlen, newunit

  implicit none

  private

  integer, parameter :: nTimers = 8

  !> Handles for timer objects to measure the time for some code parts
  integer, public :: timer_handle_loadconfig
  integer, public :: timer_handle_proto
  integer, public :: timer_handle_flooding
  integer, public :: timer_handle_inHeritDR
  integer, public :: timer_handle_refineLeaf
  integer, public :: timer_handle_smoothLeaf
  integer, public :: timer_handle_proto2treelm
  integer, public :: timer_handle_dumping_mesh

  integer :: first_timer

  public :: sdr_addTimers
  public :: sdr_dumptimers


contains


  !> Setup timers to assess the runtime of various parts of Seeder
  subroutine sdr_addTimers()
    !--------------------------------------------------------------------------
    !--------------------------------------------------------------------------

    ! Create seeder timers
    call tem_addTimer( timerHandle = timer_handle_loadconfig,  &
      &                timerName   = 'loadconfig'              )
    first_timer = timer_handle_loadconfig

    call tem_addTimer( timerHandle = timer_handle_proto,       &
      &                timerName   = 'proto'                   )
    call tem_addTimer( timerHandle = timer_handle_flooding,    &
      &                timerName   = 'flooding'                )
    call tem_addTimer( timerHandle = timer_handle_inHeritDR,   &
      &                timerName   = 'inHeritDR'               )
    call tem_addTimer( timerHandle = timer_handle_refineLeaf,  &
      &                timerName   = 'refineLeaf'              )
    call tem_addTimer( timerHandle = timer_handle_smoothLeaf,  &
      &                timerName   = 'smoothLeaf'              )
    call tem_addTimer( timerHandle = timer_handle_proto2treelm,&
      &                timerName   = 'proto2treelm'            )
    call tem_addTimer( timerHandle = timer_handle_dumping_mesh,&
      &                timerName   = 'dump_mesh'               )
  end subroutine sdr_addTimers


  !> Performance results are written to a file for statistical review
  subroutine sdr_dumptimers(general, nFluids, nBnds)
    !--------------------------------------------------------------------------
    !> Parameters of the current simulation
    type(tem_general_type), intent(in) :: general
    !> Number of fluid elements in dumped mesh
    integer, intent(in) :: nFluids
    !> Number of fluid elements which has boundary in dumped mesh
    integer, intent(in) :: nBnds
    !--------------------------------------------------------------------------
    logical                        :: file_exists
    character(len=pathLen)         :: filename
    integer                        :: fileunit, iTimer
    real(kind=rk),allocatable      :: timerVal(:)
    character(len=40),allocatable  :: timerLabel(:)
    character(len=PathLen)         :: header
    character(len=2**8)            :: output
    integer                        :: thandle
    real(kind=rk)                  :: MFEPS
    real(kind=rk)                  :: totaltime
    !--------------------------------------------------------------------------
    ! allocate the arrays for the timervalue and label
    allocate( timerVal( nTimers ))
    allocate( timerLabel( nTimers ))

    ! compute MFEPS (million fluids elements per second)
    totaltime = tem_getMaxtimerval( timerhandle = 1,                &
      &                             comm        = general%proc%comm )
    MFEPS = nFluids/(totaltime*10**6)

    ! writes the first part of the header
    write(header,'(a1,1x,a12,1x,a20,1x,a12,1x,a12,1x,a8,1x,a8)') &
            & '#', 'Revision', &
            & 'Casename',      &
            & 'nFluids',       &
            & 'nBnds',         &
            & 'nProcs',        &   ! The number of proc (i.e. MPI ranks)
            & 'MFEPS'              ! million fluids elements per second

    ! run over all timers and get the value and the label of the timer
    do iTimer = 1, nTimers
      thandle = first_timer+iTimer-1
      timerVal(iTimer) = tem_getMaxtimerval(timerhandle = thandle,          &
        &                                   comm        = general%proc%comm )
      timerLabel(iTimer) = trim(tem_getTimerName(timerHandle = thandle))
      write(header,'(a,a12,a1)') trim(header), trim(timerLabel(iTimer)), '|'
    enddo


    ! writes the first part of the output
    write(output,'(1x,a13,1x,a20,1x,i12,1x,i12,1x,i8,1x,EN12.3)') &
      &  trim(soi_solver_revision),       &
      &  trim(general%solver%simName),    &
      &  nFluids,                         &
      &  nBnds,                           &
      &  general%proc%comm_size,          &
      &  MFEPS

    ! loop to write timer value in output
    do iTimer = 1, nTimers
      write(output,'(a,1x,EN12.3)') trim(output), timerVal( iTimer )
    enddo

    !> write file
    filename = trim(general%timingFile)

    inquire(file=filename, exist=file_exists)

    fileunit = newunit()
    open(unit=fileunit, file=trim(filename), position='append')

    if (.not. file_exists ) then
       write(fileunit,'(a)') trim(header)
    end if

    write(fileunit,'(a)') trim(output)
    close(fileunit)

  end subroutine sdr_dumptimers

end module sdr_timer_module
