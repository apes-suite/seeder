! Copyright (c) 2011 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2011-2014 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2011-2012, 2014 Harald Klimach <harald.klimach@uni-siegen.de>
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
! ******************************************************************************!
!> Some auxilary functionalities.
module sdr_aux_module

  ! treelm modules
  use tem_aux_module,       only: tem_print_execInfo
  use tem_logging_module,   only: logunit
  use tem_solveHead_module, only: tem_solveHead_type

  implicit none

  private

  public :: sdr_init_global

contains

  ! ****************************************************************************!
  !> Prominently let the user now, what he actually is running right now.
  !!
  !! Also set the solvername and version number in the solveHead.
  subroutine sdr_init_global(solveHead)
    !> contains solver header information
    type( tem_solveHead_type ), intent(in) :: solveHead

    write(logunit(1),*) "                                             "
    write(logunit(1),*) "                              _              "
    write(logunit(1),*) "            ___  ___  ___  __| | ___ _ __    "
    write(logunit(1),*) "           / __|/ _ \/ _ \/ _` |/ _ \ '__|   "
    write(logunit(1),*) "           \__ \  __/  __/ (_| |  __/ |      "
    write(logunit(1),*) "           |___/\___|\___|\__,_|\___|_"      &
                        // trim(solveHead%version)
    write(logunit(1),*) "                                             "
    write(logunit(1),*) &
      &  " (C) 2012 German Research School for Simulation Sciences"
    write(logunit(1),*) " (C) 2013 University of Siegen               "
    write(logunit(1),*) "                                             "
    ! Write the information about the executable, gathered at build time to
    ! the screen.
    call tem_print_execInfo()
    write(logunit(1),*) "                                             "
    write(logunit(1),*) "                                             "

  end subroutine sdr_init_global
  ! ****************************************************************************!


end module sdr_aux_module
! ******************************************************************************!
