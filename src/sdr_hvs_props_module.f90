! Copyright (c) 2015-2016 Harald Klimach <harald.klimach@uni-siegen.de>
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
!> This module describes the properties in a treelmesh to allow their output.
module sdr_hvs_props_module
  use, intrinsic :: iso_c_binding, only: c_loc, c_f_pointer

  use env_module,               only: rk
  use treelmesh_module,         only: treelmesh_type
  use tem_bc_prop_module,       only: tem_bc_prop_type, init_tem_bc_prop
  use tem_color_prop_module,    only: tem_color_prop_type, &
    &                                 tem_color_prop_load, &
    &                                 colors_per_char
  use tem_comm_env_module,      only: tem_comm_env_type
  use tem_logging_module,       only: logunit
  use tem_time_module,          only: tem_time_type
  use tem_varsys_module,        only: tem_varSys_proc_element,             &
    &                                 tem_varSys_proc_point,               &
    &                                 tem_varSys_proc_getParams,           &
    &                                 tem_varSys_proc_setupIndices,        &
    &                                 tem_varSys_proc_getValOfIndex,       &
    &                                 tem_varsys_append_dervar,            &
    &                                 tem_varSys_type, tem_varSys_op_type, &
    &                                 tem_varSys_getParams_dummy

  use ply_dof_module,           only: Q_space
  use ply_subresolution_module, only: ply_subresolution_type,   &
    &                                 ply_subresolution_load,   &
    &                                 ply_subres_import_color,  &
    &                                 ply_subres_get_elemcolor, &
    &                                 ply_subres_colvar_type

  implicit none

  private

  !> Auxilary data type to allow an array of pointers to colvar objects.
  type sdr_hvs_colvar_ptr_type
    type(ply_subres_colvar_type), pointer :: p => NULL()
    procedure(tem_varSys_proc_element), nopass, pointer :: get_element
  end type sdr_hvs_colvar_ptr_type


  !> Datatype to collect the various properties of the mesh.
  type sdr_hvs_props_type
    !> Property describing the boundary conditions in the mesh.
    type(tem_bc_prop_type), pointer :: bc => NULL()

    !> Property describing the coloring of the mesh.
    type(tem_color_prop_type), pointer :: color => NULL()

    !> Additional color distribution information via subresolution.
    type(ply_subresolution_type), pointer :: subres => NULL()

    !> Description of color variables.
    type(sdr_hvs_colVar_ptr_type), allocatable :: colorVar(:)

    !> Position in the variables for each color.
    !!
    !! If a color was not added as variable, the entry will be negative.
    !! This is useful to iterate through all colors, that were actually added
    !! as variables.
    integer, allocatable :: colorvarpos(:)
  end type sdr_hvs_props_type

  public :: sdr_hvs_props_type
  public :: sdr_hvs_props_load
  public :: sdr_hvs_props_import_dofs


contains


  !----------------------------------------------------------------------------!
  !> Load the properties from the mesh.
  subroutine sdr_hvs_props_load(me, varsys, mesh, proc)
    !> Properties to load.
    type(sdr_hvs_props_type), intent(out) :: me

    !> A variable system, to which the colors should be appended.
    type(tem_varsys_type), intent(inout) :: varsys

    !> Mesh to load the properties for.
    type(treelmesh_type), intent(in) :: mesh

    !> Process description for MPI-IO
    type(tem_comm_env_type), intent(in) :: proc
    !----------------------------------------------------------------------!
    integer :: iCol
    logical :: wasAdded
    procedure(tem_varSys_proc_element), pointer :: get_element => NULL()
    procedure(tem_varSys_proc_point), pointer :: get_point => NULL()
    procedure(tem_varSys_proc_getParams), pointer :: get_params => NULL()
    procedure(tem_varSys_proc_setupIndices), pointer :: setup_indices => NULL()
    procedure(tem_varSys_proc_getValOfIndex), pointer :: get_valOfIndex &
      &                                                  => NULL()
    !----------------------------------------------------------------------!

    write(logunit(1),'(a)') 'Loading mesh properties...'

    call sdr_hvs_props_clean(me)

    allocate(me%bc)
    allocate(me%color)
    allocate(me%subres)

    write(logunit(5),'(a)') 'Init tem BC prop'
    call init_tem_bc_prop( tree   = mesh,      &
      &                    myPart = proc%rank, &
      &                    comm   = proc%comm, &
      &                    bc     = me%bc      )

    write(logunit(5),'(a)') 'Load color prop'
    call tem_color_prop_load( me     = me%color,  &
      &                       tree   = mesh,      &
      &                       myPart = proc%rank, &
      &                       comm   = proc%comm  )

    write(logunit(5),'(a)') 'Load subres prop'
    call ply_subresolution_load( me       = me%subres, &
      &                          tree     = mesh,      &
      &                          proc     = proc,      &
      &                          coloring = me%color   )

    get_element => ply_subres_get_elemcolor
    get_params => tem_varSys_getparams_dummy
    nullify(get_point)
    nullify(setup_indices)

    allocate(me%colorvar(me%color%nColors))
    allocate(me%colorvarpos(me%color%nColors))

    write(logunit(1),'(a,i0,a)') 'Appending ', me%color%nColors, &
      &                          ' to the varsys:'
    do iCol=1,me%color%nColors
      allocate(me%colorvar(iCol)%p)
      me%colorvar(iCol)%p%color => me%color
      me%colorvar(iCol)%p%subres => me%subres
      me%colorvar(iCol)%p%colpos = iCol
      call tem_varSys_append_derVar(                     &
        &    me             = varsys,                     &
        &    varname        = me%color%color_label(iCol), &
        &    nComponents    = 1,                          &
        &    method_data    = c_loc(me%colorvar(iCol)%p), &
        &    get_point      = get_point,                  &
        &    get_element    = get_element,                &
        &    get_params     = get_params,                 &
        &    setup_indices  = setup_indices,              &
        &    get_valofindex = get_valofindex,             &
        &    wasAdded       = wasAdded,                   &
        &    pos            = me%colorvarpos(iCol)        )

      if (.not. wasAdded) then
        ! The color was not added, due to an already existing variable with the
        ! same name, issue a warning and mark this color to skip when iterating
        ! through all colors.
        me%colorvarpos(iCol) = -1
        write(logunit(1),*) 'WARNING: Color ', &
          &                 trim(me%color%color_label(iCol)), &
          &                 ' was NOT appended!'
        write(logunit(1),*) 'Colors that have the same name as meshinfos'
        write(logunit(1),*) 'can not be considered!'
        write(logunit(1),*) ''
        write(logunit(1),*) 'Please be aware, that this name will not refer'
        write(logunit(1),*) 'to the color in the final visualization file!'
      else
        write(logunit(1),*) 'Appended Color ', &
          &                 trim(me%color%color_label(iCol)), &
          &                 ' to the available variables.'
      end if
    end do

    write(logunit(1),'(a)') '... done loading mesh properties.'

  end subroutine sdr_hvs_props_load
  !----------------------------------------------------------------------------!
  !----------------------------------------------------------------------------!


  !----------------------------------------------------------------------------!
  !> Import subresolved color data.
  subroutine sdr_hvs_props_import_dofs(me, mesh, proc, maxdegree, ndims)
    !> Properties to get the subresolution dofs for.
    type(sdr_hvs_props_type), intent(inout) :: me

    !> Mesh definition, these properties belong to.
    type(treelmesh_type), intent(in) :: mesh

    !> Communicating environment.
    type(tem_comm_env_type), intent(in) :: proc

    !> Maximal polynomial degree to import.
    !!
    !! Set to 0, if only the integral mean value is to be loaded.
    integer, intent(in) :: maxdegree

    !> Dimensions of the polynomials after importing.
    integer, intent(in) :: ndims
    !----------------------------------------------------------------------!
    integer :: iColor, cpos
    !----------------------------------------------------------------------!

    do iColor=1,me%color%nColors

      cpos = me%colorvarpos(iColor)

      if (cpos > 0) then
        ! Skip colors that where not added to the list of variables. And only
        ! process those with a valid variable position.

        ! Import the subresolution data for each color and store it in the
        ! subresdat of each variable method data.
        ! We always use Q_space here.
        call ply_subres_import_color( me            = me%subres,          &
          &                           tree          = mesh,               &
          &                           iColor        = iColor,             &
          &                           coloring      = me%color,           &
          &                           target_degree = maxdegree,          &
          &                           target_space  = Q_space,            &
          &                           target_dim    = ndims,              &
          &                           subresdat     = me%colorvar(iColor) &
          &                                             %p%subresdat      )
        me%colorvar(iColor)%p%nsubdofs = (maxdegree+1)**nDims
      end if

    end do

  end subroutine sdr_hvs_props_import_dofs
  !----------------------------------------------------------------------------!
  !----------------------------------------------------------------------------!


  !----------------------------------------------------------------------------!
  !> Clean up the properties data type.
  !!
  !! Takes care of the pointers, nullifies and deallocates them properly.
  subroutine sdr_hvs_props_clean(me)
    type(sdr_hvs_props_type), intent(out) :: me
    !----------------------------------------------------------------------!
    integer :: iCol
    !----------------------------------------------------------------------!

    if (allocated(me%colorvar)) then
      do iCol=lbound(me%colorvar,1),ubound(me%colorvar,1)
        if (associated(me%colorvar(iCol)%p)) then
          nullify(me%colorvar(iCol)%p%color)
          nullify(me%colorvar(iCol)%p%subres)
          if ( allocated(me%colorvar(iCol)%p%subresdat) ) then
            deallocate(me%colorvar(iCol)%p%subresdat)
          end if
          deallocate(me%colorvar(iCol)%p)
        end if
      end do
      deallocate(me%colorvar)
    end if

    if (allocated(me%colorvarpos)) then
      deallocate(me%colorvarpos)
    end if

    if (associated(me%bc)) deallocate(me%bc)
    if (associated(me%color)) deallocate(me%color)
    if (associated(me%subres)) deallocate(me%subres)

  end subroutine sdr_hvs_props_clean

end module sdr_hvs_props_module
