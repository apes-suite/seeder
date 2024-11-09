! Copyright (c) 2014 Harald Klimach <harald.klimach@uni-siegen.de>
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
module sdr_subres_fills_module
  use env_module,             only: rk, labelLen
  use flu_binding,            only: flu_State
  use aotus_module,           only: aot_get_val, flu_State, &
    &                               aoterr_Fatal
  use aot_table_module,       only: aot_table_open, aot_table_close, &
    &                               aot_table_length

  use tem_aux_module,         only: tem_abort
  use tem_logging_module,     only: logunit
  use tem_tools_module,       only: upper_to_lower
  use tem_dyn_array_module,   only: dyn_labelArray_type, init, append
  use tem_grow_array_module,  only: grw_realArray_type, init, append

  implicit none

  !> Definition of values to use for colors in polynomial representations.
  type sdr_subres_fills_type
    !> Name of the color, to which the color values apply.
    type(dyn_labelArray_type) :: label

    !> Actual value to use, where the color is present, defaults to 1.
    !!@todo HK: replace by spatial function
    type(grw_realArray_type) :: fill

    !> Actual value to use, where the color is not present, defaults to 0.
    type(grw_realArray_type) :: void
  end type sdr_subres_fills_type

  public :: sdr_subres_fills_type
  public :: sdr_subres_fills_load
  public :: sdr_subres_fills_add


contains


  ! ****************************************************************************!
  !> Add a value definition for a color to the list of fillings.
  subroutine sdr_subres_fills_add( fills, colorname, fill_value, void_value, &
    &                              pos, wasAdded )
    ! --------------------------------------------------------------------------!
    !> Table of color values to add a color definition to.
    type(sdr_subres_fills_type), intent(inout) :: fills
    !> Name of the color to define values for.
    character(len=*), intent(in) :: colorname
    !> Value to use, where the domain is flooded by color.
    real(kind=rk), intent(in) :: fill_value
    !> Value to use, where the domain is NOT flooded by color.
    real(kind=rk), intent(in) :: void_value
    !> Position of this color in the list of color values.
    integer, intent(out) :: pos
    !> Indicator, if this value definition was added to the list.
    logical, intent(out) :: wasAdded
    ! --------------------------------------------------------------------------!
    ! --------------------------------------------------------------------------!

    call append(fills%label, upper_to_lower(trim(colorname)), &
      &         pos=pos, wasAdded = wasAdded)

    if (wasAdded) then
      call append(fills%fill, fill_value)
      call append(fills%void, void_value)
    end if

  end subroutine sdr_subres_fills_add
  ! ****************************************************************************!


  ! ****************************************************************************!
  !> Load the filling definition for subresolved colors.
  !!
  !! A filling is defined by a value for the area where the domain is flooded
  !! by the color specified in the label and a value for the remaining part
  !! of the domain without flooding.
  !!
  !! The filling information has to be provided as a list of tables in the
  !! following form within the subresolution definition:
  !! \code{.lua}
  !! values = {
  !!            { label = 'colA', fill = 1.0, void = 0.0 },
  !!            { label = 'colB', fill = 9.0, void = 0.1 }
  !!          }
  !! \endcode
  !! Defaults are: fill=1.0, void=0.0.
  subroutine sdr_subres_fills_load(fills, conf, parent)
    ! --------------------------------------------------------------------------!
    !> Value definitions for individual colors.
    type(sdr_subres_fills_type), intent(out) :: fills

    !> Handle for the Lua script.
    type(flu_State) :: conf

    !> Parent table, within which to open the fills table.
    integer, intent(in) :: parent
    ! --------------------------------------------------------------------------!
    integer :: thandle
    integer :: nValues
    integer :: iValue
    ! --------------------------------------------------------------------------!

    call aot_table_open(L = conf, parent = parent, thandle = thandle, &
      &                 key = 'values')

    if (thandle > 0) then
      nValues = aot_table_length(L = conf, thandle = thandle)
      call init(fills%label, length = nValues)
      call init(fills%fill, length = nValues)
      call init(fills%void, length = nValues)
      do iValue=1,nValues
        call load_single_fill( fills  = fills,   &
          &                    conf   = conf,    &
          &                    parent = thandle, &
          &                    pos    = iValue   )
      end do
    end if

  end subroutine sdr_subres_fills_load
  ! ****************************************************************************!


  ! ****************************************************************************!
  !> Load the value fill definition for a single color.
  !!
  !!
  !! The filling information has to be provided as a table with the label of
  !! the color, the fill value to use, where the color is present and the void
  !! value for the rest of the domain:
  !! \code{.lua}
  !!            { label = 'colA', fill = 1.0, void = 0.0 },
  !! \endcode
  !!
  !! Defaults are: fill=1.0, void=0.0.
  !! Note that the keys can be left out, and assignment by position is possible.
  subroutine load_single_fill(fills, conf, parent, pos)
    ! --------------------------------------------------------------------------!
    !> Value definitions for individual colors.
    type(sdr_subres_fills_type), intent(inout) :: fills

    !> Handle for the Lua script.
    type(flu_State) :: conf

    !> Parent table, within which to open the fills table.
    integer, intent(in) :: parent

    !> Position in the table of fills to load.
    integer, intent(in) :: pos
    ! --------------------------------------------------------------------------!
    integer :: thandle
    character(len=labelLen) :: cLabel
    real(kind=rk) :: fill_value
    real(kind=rk) :: void_value
    integer :: iError
    integer :: colpos
    logical :: wasAdded
    ! --------------------------------------------------------------------------!

    call aot_table_open(L = conf, parent = parent, thandle = thandle, &
      &                 pos = pos)

    ! The label of the color for which we set the value is either in position 1
    ! or in the key label.
    call aot_get_val( L       = conf,    &
      &               thandle = thandle, &
      &               key     = 'label', &
      &               pos     = 1,       &
      &               val     = clabel,  &
      &               ErrCode = iError   )

    if (btest(iError, aoterr_Fatal)) then
      write(logunit(0),*) &
        &  'FATAL Error occured while retrieving a filling color name!'
      call tem_abort()
    end if

    ! The filling value is either found with the key fill, or as second
    ! argument.
    call aot_get_val( L       = conf,       &
      &               thandle = thandle,    &
      &               key     = 'fill',     &
      &               pos     = 2,          &
      &               val     = fill_value, &
      &               default = 1.0_rk,     &
      &               ErrCode = iError      )

    if (btest(iError, aoterr_Fatal)) then
      write(logunit(0),*) &
        &  'FATAL Error occured while retrieving a filling value for color ' &
        &  // trim(clabel)
      call tem_abort()
    end if


    ! The void value is either found with the key void, or as third argument.
    call aot_get_val( L       = conf,       &
      &               thandle = thandle,    &
      &               key     = 'void',     &
      &               pos     = 3,          &
      &               val     = void_value, &
      &               default = 0.0_rk,     &
      &               ErrCode = iError      )

    if (btest(iError, aoterr_Fatal)) then
      write(logunit(0),*) &
        &  'FATAL Error occured while retrieving a void value for color ' &
        &  // trim(clabel)
      call tem_abort()
    end if

    call aot_table_close(L = conf, thandle = thandle)

    ! Now add the color values to the list.
    call sdr_subres_fills_add( fills      = fills,      &
      &                        colorname  = cLabel,     &
      &                        fill_value = fill_value, &
      &                        void_value = void_value, &
      &                        pos        = colpos,     &
      &                        wasAdded   = wasAdded    )

    if (.not. wasAdded) then
      write(logunit(1),*) 'WARNING: subresolution values for color ', cLabel
      write(logunit(1),*) 'already in the list!'
      write(logunit(1),*) 'Keep in mind that color names are case' &
        &                 // ' insensitive.'
    end if

  end subroutine load_single_fill
  ! ****************************************************************************!

end module sdr_subres_fills_module
