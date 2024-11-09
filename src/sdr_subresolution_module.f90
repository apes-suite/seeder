! Copyright (c) 2014-2015 Harald Klimach <harald.klimach@uni-siegen.de>
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
!> This module describes settings that are required for the resolution of
!! boundaries on a subelement level.
!!
!! These boundaries are represented by polynomials, and a nodal sampling in
!! elements intersected by the boundary needs to be done accordingly.
!!
!! For subresolution to work, there needs to be a subresolution table
!! defined, indicating at least the polynomial degree to use for polynomial
!! representations.
!! Only boundaries with a color (not the 'none' color), can be subresolved.
!! For each color a separate file will be generated to contain the
!! polynomial information on all the subresolved elements of this color.
!! Note, that polynomial settings are global, all boundaries will be resolved
!! with the same polynomial definition.
!!
!! Nodal values are obtained, by refining elements further down and probing
!! integration points for their in or out status. This detaches the
!! actual polynomial somewhat from the mesh resolution, and the number of
!! levels to resolve beyond the target element can be stated independently
!! of the polynomial degree.
!! However, the subresolution levels should be set sufficiently high to
!! provide a good resolution for the integration points.
module sdr_subresolution_module
  use env_module,                  only: rk, labelLen
  use flu_binding,                 only: flu_State
  use aotus_module,                only: aot_get_val, flu_State, &
    &                                    aoterr_Fatal
  use aot_table_module,            only: aot_table_open, aot_table_close, &
    &                                    aot_table_length

  use ply_dof_module,              only: P_Space, Q_Space
  use ply_dynarray_project_module, only: ply_prj_init_type, ply_prj_init_define
  use ply_prj_header_module,       only: ply_prj_header_type, &
    &                                    ply_prj_header_load
  use tem_aux_module,              only: tem_abort
  use tem_logging_module,          only: logunit
  use tem_tools_module,            only: upper_to_lower
  use tem_dyn_array_module,        only: dyn_labelArray_type, init, append
  use sdr_subres_fills_module,     only: sdr_subres_fills_type, &
    &                                    sdr_subres_fills_load, &
    &                                    sdr_subres_fills_add

  implicit none

  !> Settings for the resolution below mesh elements.
  type sdr_subresolution_type
    !> Degree for the polynomials to use in the subelement resolution.
    integer :: polydegree

    !> Identification of the polynomial space to use.
    character :: polyspace

    !> Definition of conversions between modal and nodal values.
    type(ply_prj_init_type) :: projection

    !> Header definition of the projection header.
    !!
    !! Stored here for later output in the subresolution file.
    type(ply_prj_header_type) :: project_header

    !> Number of levels to use for the resolution of boundaries within elements.
    !!
    !! By default this value is automatically determined by the polydegree and
    !! the point set that is to be used. But that can be explicitly overwritten
    !! by the user with an arbitrary value here.
    !! If there is no subresolution to be done, this value will be 0.
    integer :: nLevels

    !> Dynamic array to store the color names, for which the default boundary
    !! resolution should reach subelements.
    !!
    !! All boundaries with this color will by default have subelement
    !! resolution, but this default can still be overwritten by the attribute
    !! setting of the boundary.
    type(dyn_labelArray_type) :: default_colors

    !> Definition of the values to use for color filling and color voids.
    type(sdr_subres_fills_type) :: color_values

    !> Position of value definitions for each color in the mesh.
    integer, allocatable :: value_pos(:)
  end type sdr_subresolution_type

  public :: sdr_subresolution_type
  public :: sdr_subresolution_load
  public :: sdr_subresolution_encolor


contains


  !> Loading the settings for the subresolution in the mesh.
  subroutine sdr_subresolution_load(subres, conf, parent)
    !> The subresolution data structure to fill.
    type(sdr_subresolution_type), intent(out) :: subres

    !> Handle to the Lua configuration script.
    type(flu_State) :: conf

    !> Handle to a possible parent table.
    integer, optional, intent(in) :: parent

    integer :: basisType
    integer :: thandle
    integer :: phandle
    integer :: iError
    integer :: deflev
    character(len=labelLen) :: polyspace
    real(kind=rk) :: ofact

    ! Default: no subelement resolution:
    subres%polydegree = 0
    subres%nLevels = 0

    call aot_table_open(L = conf, parent = parent, thandle = thandle, &
      &                 key = 'subresolution')

    if (thandle > 0) then
      call aot_get_val( L       = conf,              &
        &               thandle = thandle,           &
        &               key     = 'polydegree',      &
        &               val     = subres%polydegree, &
        &               ErrCode = iError             )

      if (btest(iError, aoterr_Fatal)) then
        write(logUnit(0),*) &
          &  'FATAL Error occured, while retrieving subresolution polydegree'
        call tem_abort()
      end if

      ! All other settings are only useful, if the polydegree is actually
      ! larger than 0!
      if (subres%polydegree > 0) then
        call aot_get_val( L       = conf,        &
          &               thandle = thandle,     &
          &               key     = 'polyspace', &
          &               val     = polyspace,   &
          &               default = 'q',         &
          &               ErrCode = iError       )

        select case(upper_to_lower(trim(polyspace)))
        case('q')
          basisType = Q_space

        case('p')
          basisType = P_space

        case default
          write(logunit(1),*) 'ERROR in subresolution loading!'
          write(logunit(1),*) 'Unknown polyspace ', trim(polyspace)
          write(logUnit(1),*) 'Supported are:'
          write(logUnit(1),*) '* Q (quadratic with i,j,k <= maxDegree)'
          write(logUnit(1),*) '* P (with i+j+k <= maxDegree)'
          write(logUnit(1),*) 'Stopping....'
          call tem_abort()

        end select

        subres%polyspace = upper_to_lower(trim(polyspace))

        call aot_table_open(L = conf, parent = thandle, thandle = phandle, &
          &                 key = 'projection')
        call ply_prj_header_load( me     = subres%project_header, &
          &                       conf   = conf,                  &
          &                       parent = phandle                )

        if (trim(subres%project_header%kind) == 'fpt') then
          ofact = subres%project_header%fpt_header%factor
        else
          ofact = subres%project_header%l2p_header%factor
        end if

        call aot_table_close(L = conf, thandle = phandle)

        call ply_prj_init_define( me            = subres%projection,     &
          &                       header        = subres%project_header, &
          &                       maxPolyDegree = subres%polydegree,     &
          &                       basisType     = basisType              )

        deflev = ceiling(2*log(real((subres%polydegree+1)*ofact))/log(2.0))
        write(logunit(2),*) 'Reading the subresolving level to use.'
        write(logunit(2),*) 'Default: ', deflev
        call aot_get_val(L = conf, thandle = thandle, key = 'levels', &
          &              val = subres%nLevels, default = deflev, &
          &              ErrCode = iError )

        if (btest(iError, aoterr_Fatal)) then
          write(logUnit(0),*) &
            &  'FATAL Error occured, while retrieving subresolution levels'
          call tem_abort()
        end if

        if (subres%nLevels < 1) then
          write(logUnit(0),*) 'If a subelement resolution by polynomials should'
          write(logUnit(0),*) 'be done, levels needs to be at least 1!'
          write(logUnit(0),*) 'But levels = ', subres%nLevels
          write(logUnit(0),*) 'FATAL Error! Check your configuration.'
          call tem_abort()
        end if

        call subres_load_color_defaults( label  = subres%default_colors, &
          &                              conf   = conf,                  &
          &                              parent = thandle                )

      end if

      ! Fills might even be defined, if there is no polynomial for
      ! subresolution. It still can be used to define color values.
      call sdr_subres_fills_load( fills  = subres%color_values, &
        &                         conf   = conf,                &
        &                         parent = thandle              )

    else

      ! Initialize the subres fills, also if there is no subresolution table.
      ! This is a little counter-intuitive, but we need to provide color values
      ! in any case, thus the arrays need to be properly initialized here.
      call sdr_subres_fills_load( fills  = subres%color_values, &
        &                         conf   = conf,                &
        &                         parent = parent               )

    end if

    call aot_table_close(L = conf, thandle = thandle)

  end subroutine sdr_subresolution_load



  !> Internal subroutine to load the list of colors, which by default should
  !! apply subelement resolution to its boundaries.
  !!
  !! This list is provided in the subresolution table by:
  !! \code{.lua}
  !! default_colors = { 'colA', 'colB' }
  !! \endcode
  !! This table is read by this routine.
  subroutine subres_load_color_defaults(label, conf, parent)
    !> Labels of the colors, for which a default setting is given.
    type(dyn_labelArray_type), intent(out) :: label

    !> Handle for the Lua script.
    type(flu_State) :: conf

    !> Parent table, within which to open the color_default table.
    integer, intent(in) :: parent

    character(len=labelLen) :: cLabel
    integer :: thandle
    integer :: nLabels
    integer :: iLabel
    integer :: iError
    integer :: cPos
    logical :: wasAdded

    call aot_table_open(L = conf, parent = parent, thandle = thandle, &
      &                 key = 'default_colors')

    if (thandle > 0) then
      nLabels = aot_table_length(L = conf, thandle = thandle)
      call init(label, length = nLabels)

      do iLabel=1,nLabels
        call aot_get_val(L       = conf,    &
          &              thandle = thandle, &
          &              pos     = iLabel,  &
          &              val     = clabel,  &
          &              ErrCode = iError   )

        if (btest(iError, aoterr_Fatal)) then
          write(logunit(0),*) &
            &  'FATAL Error occured while retrieving a subresolved default' &
            &  // ' color name!'
          call tem_abort()
        end if

        call append(label, upper_to_lower(trim(clabel)), &
          &         pos=cPos, wasAdded = wasAdded)

        if (.not. wasAdded) then
          write(logunit(1),*) 'WARNING: subresolution default color ', cLabel
          write(logunit(1),*) 'already in the list!'
          write(logunit(1),*) 'Keep in mind that color names are case' &
            &                 // 'insensitive'
        end if
      end do
    end if

    call aot_table_close(L = conf, thandle = thandle)

  end subroutine subres_load_color_defaults


  !> Find the value definitions for all unique colors.
  subroutine sdr_subresolution_encolor(me, colors)
    type(sdr_subresolution_type), intent(inout) :: me
    type(dyn_labelArray_type), intent(in) :: colors

    integer :: iColor
    character(len=labelLen) :: cLabel
    logical :: wasAdded

    allocate(me%value_pos(colors%nVals))
    me%value_pos = -1

    do iColor=1,colors%nVals
      cLabel = colors%val(iColor)
      if (trim(cLabel) /= 'none') then
        write(logunit(4), *) 'Reading filling definition for color ' &
          &                  //trim(cLabel)
        call sdr_subres_fills_add( fills      = me%color_values,      &
          &                        colorname  = trim(cLabel),         &
          &                        fill_value = 1.0_rk,               &
          &                        void_value = 0.0_rk,               &
          &                        pos        = me%value_pos(iColor), &
          &                        wasAdded   = wasAdded              )
        if (wasAdded) then
          write(logunit(1), *) 'NOTE: No values specified for color ', &
            &                  trim(cLabel)
          write(logunit(1), *) 'Using default fill=1, void=0 !!!'
        end if
        write(logunit(2), *) 'Color values for ' // trim(cLabel) &
          &                  // ': fill=', &
          &                  me%color_values%fill%val(me%value_pos(iColor)), &
          &                  '; void=', &
          &                  me%color_values%void%val(me%value_pos(iColor))
      end if
    end do

  end subroutine sdr_subresolution_encolor

end module sdr_subresolution_module
