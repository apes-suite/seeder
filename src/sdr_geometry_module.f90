! Copyright (c) 2012-2015, 2017 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2012-2015, 2019 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2012, 2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
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
!> This module provides the description of the complete geometry to mesh.
module sdr_geometry_module
  use env_module,                only: rk, labelLen
  use tem_aux_module,            only: tem_abort
  use tem_tools_module,          only: upper_to_lower, tem_horizontalSpacer
  use tem_logging_module,        only: logunit
  use tem_dyn_array_module,      only: dyn_labelArray_type, PositionOfVal
  use tem_grow_array_module,     only: grw_intArray_type
  use tem_cube_module,           only: tem_cube_type, tem_load_cube
  use tem_triangle_module,       only: grw_triangleArray_type, init, truncate, &
    &                                  tem_triangleCubeOverlap,                &
    &                                  tem_load_triangle
  use tem_point_module,          only: grw_pointArray_type, init, truncate,    &
    &                                  tem_pointCubeOverlap
  use tem_line_module,           only: grw_lineArray_type, init, truncate,     &
    &                                  tem_lineCubeOverlap
  use tem_box_module,            only: grw_boxArray_type, init, truncate,      &
    &                                  tem_boxCubeOverlap
  use tem_sphere_module,         only: tem_load_sphere, tem_sphereCubeOverlap, &
    &                                  grw_sphereArray_type, init, truncate,   &
    &                                  append, tem_sphere_type
  use tem_ellipsoid_module,      only: tem_load_ellipsoid,                     &
    &                                  tem_ellipsoidCubeOverlap,               &
    &                                  grw_ellipsoidArray_type, init, truncate,&
    &                                  append, tem_ellipsoid_type
  use tem_cylinder_module,       only: tem_load_cylinder,                      &
    &                                  tem_cylinderCubeOverlap,                &
    &                                  grw_cylinderArray_type, init, truncate
  use tem_stl_module,            only: tem_load_stl
  use tem_transformation_module, only: tem_transformation_type,                &
    &                                  tem_load_transformation


  use sdr_spatialObj_module,     only: grw_spatialObjArray_type,               &
    &                                  sdr_spatialObj_type,                    &
    &                                  point, line, triangle, box, sphere, ellipsoid,&
    &                                  spacerInterwoven, periodicPlane,        &
    &                                  cylinder, init, append, truncate
  use sdr_periodic_module,       only: grw_periPlaneArray_type,                &
    &                                  init, truncate,                         &
    &                                  sdr_load_periodic,                      &
    &                                  sdr_periodicPlaneCubeOverlap
  use sdr_attribute_module,      only: sdr_load_attribute, sdr_init_attribute, &
    &                                  sdr_attrList_type, init, truncate,      &
    &                                  append, sdr_attribute_type,             &
    &                                  sdr_Refinement_object,                  &
    &                                  sdr_Boundary_object,                    &
    &                                  sdr_Seed_object,                        &
    &                                  sdr_identify_bc_colors,                 &
    &                                  sdr_identify_inv_colors,                &
    &                                  sdr_isPeriodicDefined
  use sdr_spacer_module,         only: sdr_load_spacer,                        &
    &                                  sdr_spacerInterwovenCubeOverlap,        &
    &                                  grw_spacerInterwovenArray_type,         &
    &                                  init, truncate
  use sdr_canonicalND_module,    only: sdr_load_canonicalND  
  use sdr_triangle_module,       only: sdr_load_triangle
  use sdr_stl_module,            only: sdr_load_stl
  use sdr_sphere_module,         only: sdr_load_sphere
  use sdr_ellipsoid_module,      only: sdr_load_ellipsoid
  use sdr_cylinder_module,       only: sdr_load_cylinder

  use flu_binding,               only: flu_State
  use aotus_module,              only: aot_get_val, flu_State,                 &
    &                                  aoterr_Fatal, aoterr_NonExistent,       &
    &                                  aoterr_WrongType
  use aot_table_module,          only: aot_table_open, aot_table_close,        &
    &                                  aot_table_length, aot_table_top,        &
    &                                  aot_table_first

  implicit none

  public :: sdr_geometry_type
  public :: sdr_load_geometry
  public :: is_intersecting
  public :: sdr_append_distanceRefineObject


  !> Geometric description of the space to mesh.
  type sdr_geometry_type
    !> The universe given as a bounding cube.
    type(tem_cube_type) :: universe

    !> growing array of all points
    type(grw_pointArray_type) :: point

    !> growing array of all lines
    type(grw_lineArray_type) :: line

    !> growing array of all triangles
    type(grw_triangleArray_type) :: triangle

    !> growing array of all boxes
    type(grw_boxArray_type) :: box

    !> growing array of periodic plane objects
    type(grw_periPlaneArray_type) :: periPlane

    !> growing array of spheres
    type(grw_sphereArray_type) :: sphere

    !> growing array of spheres
    type(grw_ellipsoidArray_type) :: ellipsoid

    !> growing array of cylinders
    type(grw_cylinderArray_type) :: cylinder

    !> growing array of spacer
    type(grw_spacerInterwovenArray_type) :: spacerInterwoven

    !> Growing array of geometrical objects.
    type(grw_spatialObjArray_type) :: spatialObj

    !> Dynamic Unique array of attributes connected to the objects
    !! (seed, boundary or refinement)
    type(sdr_attrList_type) :: attribute

    !> Smooth refinement towards boundaries.
    logical :: smoothbounds = .true.

    !> Smooth level jumps, avoid level jump > 1 in domain
    logical :: smoothLevels = .true.

    !> Number of user defined spatial objects
    !! i.e excluding spatial objects created for distance refine
    integer :: nUserObjs
  end type sdr_geometry_type

contains

  !> Test the intersection between the given cube and the object specified by
  !! obj_pos.
  !!
  !! This test obviously depends on the kind of object to intersect, thus a
  !! select case has to be used. Alternatively we might think of deploying
  !! function pointers.
  function is_intersecting(cube, geometry, obj_pos) result(intersects)
    ! --------------------------------------------------------------------------!
    type(tem_cube_type), intent(in) :: cube
    type(sdr_geometry_type), intent(in) :: geometry
    integer, intent(in) :: obj_pos
    logical :: intersects
    ! --------------------------------------------------------------------------!
    integer :: prim_pos
    integer :: prim_kind
    ! --------------------------------------------------------------------------!

    intersects = .false.

    prim_pos  = geometry%spatialObj%val(obj_pos)%primitive_position
    prim_kind = geometry%spatialObj%val(obj_pos)%geometry_primitive

    select case (prim_kind)
    case (point)
      intersects = tem_pointCubeOverlap(geometry%point%val(prim_pos), cube)

    case (line)
      intersects = tem_lineCubeOverlap(geometry%line%val(prim_pos), cube)
      
    case (triangle)
      intersects = tem_triangleCubeOverlap(geometry%triangle%val(prim_pos), &
        &                                  cube                             )

    case (box)
      intersects = tem_boxCubeOverlap(geometry%box%val(prim_pos), cube)

    case (sphere)
      intersects = tem_sphereCubeOverlap(geometry%sphere%val(prim_pos), cube)

    case (ellipsoid)
      intersects = tem_ellipsoidCubeOverlap(geometry%ellipsoid%val(prim_pos), &
        &                                   cube                              )

    case (cylinder)
      intersects = tem_cylinderCubeOverlap(geometry%cylinder%val(prim_pos), &
        &                                  cube                             )

    case (spacerInterwoven)
      intersects = sdr_spacerInterwovenCubeOverlap( &
        &                       geometry%spacerInterwoven%val(prim_pos), cube)

    case (periodicPlane)
      intersects = sdr_periodicPlaneCubeOverlap(&
        &                   geometry%periPlane%val(prim_pos), cube)

    case default
      write(logunit(1),*) 'WARNING: Unknown primitive geometry ignored'

    end select

  end function is_intersecting


  ! *****************************************************************************
  !> Routine to load spatial object defined in config file and store in 
  !! geometry type
  subroutine sdr_load_geometry( me, subres_colors, invert_colors, conf )
    ! --------------------------------------------------------------------------!
    !> contains all geometry infomation
    type( sdr_geometry_type ), intent(out) :: me

    !> All the colors, which should be default use subelement resolution for
    !! their boundaries.
    type(dyn_labelArray_type), intent(in) :: subres_colors
    type(dyn_labelArray_type), intent(in) :: invert_colors

    type(flu_State) :: conf !< lua state
    ! --------------------------------------------------------------------------!
    integer :: spa_handle, spa_subHandle
    integer :: spa_length, iSpaObj
    integer :: iColor, nColors
    integer :: iBnd, nBounds
    integer :: iError
    ! --------------------------------------------------------------------------!
    
    ! Read the bounding cube definition from the config file.
    write(logunit(1),*) 'Loading bounding cube:'
    call tem_load_cube(me = me%universe, conf = conf, key = 'bounding_cube')

    ! Smoothbounds sets, whether refinement towards boundaries should be
    ! somewhat smoothed or not.
    call aot_get_val(L=conf, key='smoothbounds', &
      &              val=me%smoothbounds, ErrCode = iError, default=.true.)
    if (me%smoothbounds) then
      write(logunit(1),*) 'Smoothed refinement towards boundaries.'
    else
      write(logunit(1),*) '*** NO Smoothed refinement towards boundaries. ***'
    end if

    ! Smoothlevels sets, whether to allow level jump > 1 in domain
    call aot_get_val(L=conf, key='smoothlevels', &
      &              val=me%smoothLevels, ErrCode = iError, default=.true.)
    if (me%smoothLevels) then
      write(logunit(1),*) 'Smoothed levels with level jump = 1.'
    else
      write(logunit(1),*) '*** NO Smoothed levels. ***'
    end if


    ! Intialize growing array of geom object.
    call init(me%spatialObj)

    ! Intialize growing arrays of geometric primitives.
    call init(me%point, length=8)
    call init(me%line, length=8)
    call init(me%triangle, length=128)
    call init(me%sphere, length=8)
    call init(me%ellipsoid, length=8)
    call init(me%cylinder, length=8)
    call init(me%spacerInterwoven, length=8)
    call init(me%periPlane, length=8)
    call init(me%box, length=8)

    ! Initialize dynamic array of attributes.
    call sdr_init_attribute(me%attribute)
    ! Load spatial objects.
    call aot_table_open(L = conf, thandle = spa_handle, key = 'spatial_object')

    ! Load spatial sub handle to check between single table or multiple table
    call aot_table_open( L = conf, parent = spa_handle,   &
      &                  thandle = spa_subhandle, pos = 1 )

    if (spa_handle .gt. 0) then
      ! Spatial Object definition exists in the configuration and is a table.

      if (spa_subHandle .eq. 0 ) then
        ! There is no second table within the opened spatial object table at the
        ! first position. Interpret the parent table as a single spatial object.
        call aot_table_close( L = conf, thandle = spa_subHandle )
        call sdr_load_spatialObject_single(me = me, conf = conf, &
          &                                subres_colors = subres_colors, &
          &                                thandle = spa_handle)
      else
        ! The first entry in the spatial object table is a table in itself, try
        ! to interpret all entries as object definitions and read multiple
        ! objects from the spatial object table.
        call aot_table_close(L = conf, thandle = spa_subHandle)
        spa_length = aot_table_length(L = conf, thandle = spa_handle)
        do iSpaObj = 1, spa_length
          call aot_table_open(L = conf, parent = spa_handle,        &
            &                 thandle = spa_subhandle, pos = iSpaObj)
          call sdr_load_spatialObject_single(me = me, conf = conf,  &
            &                                subres_colors = subres_colors, &
            &                                thandle = spa_subhandle)
          call aot_table_close(L = conf, thandle = spa_subHandle)
        end do
      endif  
    else
      write(logunit(1),*) 'WARNING: No geometrical objects are defined!'
    end if

    call aot_table_close(L = conf, thandle = spa_handle)

    ! Look up bc colors in the list of seed colors to match them correctly
    ! during flooding.
    call sdr_identify_bc_colors(me%attribute)

    ! Look up colors that are to be inverted and flag them accordingly.
    call sdr_identify_inv_colors(me%attribute, invert_colors)

    ! Truncate all arrays, to free unused memory.
    call truncate(me%spatialObj)
    call truncate(me%attribute%dynArray)
    call truncate(me%point)
    call truncate(me%line)
    call truncate(me%triangle)
    call truncate(me%sphere)
    call truncate(me%ellipsoid)
    call truncate(me%cylinder)
    call truncate(me%spacerInterwoven)
    call truncate(me%periPlane)
    call truncate(me%box)

    ! Number of user defined spatial objects to exclude spatial objects
    ! defined by distance refine
    me%nUserObjs = me%spatialObj%nVals
    write(logunit(1),*) ''
    write(logunit(1),"(A,I0)") 'Total number of spatial objs created:', &
      &                        me%nUserObjs
    write(logunit(1),"(A,I0)") 'Total number of boundaries defined:', &
      &                        me%attribute%kindpos(sdr_Boundary_object)%nVals
    nBounds = me%attribute%uni_name(sdr_Boundary_object)%nVals
    if (nBounds > 0) then
      write(logunit(1),"(A,I0,A)") 'There are ', nBounds,      &
        &                          ' unique boundaries defined:'
      do iBnd=1,nBounds
        write(logunit(1),'(A)') ' * '//trim( me%attribute                     &
          &                                    %uni_name(sdr_Boundary_object) &
          &                                    %val(iBnd)                     )
      end do
    end if
    write(logunit(1),*) ''
    nColors = me%attribute%uni_name(sdr_seed_object)%nVals
    if (nColors > 1) then
      write(logunit(1), '(a,i0,a)') 'There are ', nColors, &
        &                           ' colors in this mesh:'
      do iColor=1,nColors
        if (me%attribute%color_inverted(iColor)) then
          write(logunit(1), *) trim(me%attribute%uni_name(sdr_seed_object) &
            &                                   %val(iColor))              &
            &                  //' (inverted)'
        else
          write(logunit(1), *) trim(me%attribute%uni_name(sdr_seed_object) &
            &                                   %val(iColor))
        end if
      end do
    else
      write(logunit(1),*) 'This is a single color mesh.'
    end if
    write(logunit(1),*) ''
  
  end subroutine sdr_load_geometry
  ! *****************************************************************************


  ! *****************************************************************************
  !> Routine to load single spatial object table defined in config file
  subroutine sdr_load_spatialObject_single(me, subres_colors, conf, thandle)
    ! --------------------------------------------------------------------------!
    !> contains all geometry infomation
    type(sdr_geometry_type), intent(inout) :: me
    !> All the colors, which should be default use subelement resolution for
    !! their boundaries.
    type(dyn_labelArray_type), intent(in) :: subres_colors
    type(flu_State) :: conf !< lua state
    integer, intent(in) :: thandle !< handle for spatialObject
    ! --------------------------------------------------------------------------!
    integer :: attr_pos
    type(tem_transformation_type) :: transform
    ! --------------------------------------------------------------------------!
    call tem_horizontalSpacer(funit=logunit(1))

    ! Load the attribute of the object and append it to the list.
    call sdr_load_attribute(attrList = me%attribute, conf = conf, &
      &                     subres_colors = subres_colors, &
      &                     thandle = thandle, attr_pos = attr_pos)

    ! load transformation table if defined and pass the transformation
    ! to geom table to apply the transformation to geometry object
    call tem_load_transformation(transform = transform, conf = conf, &
      &                          thandle = thandle)

    ! Load the geometry definition of the spatial object (and connect its
    ! primitives to the current attribute).
    call sdr_load_geom_table(me = me, conf = conf, thandle = thandle, &
      &                      attr_pos = attr_pos, transform = transform)

  end subroutine sdr_load_spatialObject_single
  ! *****************************************************************************


  ! *****************************************************************************
  !> This routine loads the geometry table from the config file i.e loading
  !! different geometry kinds like canoND, cube, periodic, STL etc.
  subroutine sdr_load_geom_table(me, conf, thandle, attr_pos, transform)
    ! --------------------------------------------------------------------------!
    type(sdr_geometry_type), intent(inout) :: me
    type(flu_State) :: conf !< lua state
    integer, intent(in) :: thandle !< handle for spatialObject
    !> Position of the attribute to connect the geometric objects to.
    integer, intent(in) :: attr_pos
    type(tem_transformation_type), intent(in) :: transform
    ! --------------------------------------------------------------------------!
    integer :: geom_handle, geom_subHandle
    integer :: nTables, iTab
    ! --------------------------------------------------------------------------!

    write(logunit(1),*) 'Loading geometry '
    ! Open the geometry table.
    call aot_table_open(L = conf, parent = thandle, thandle = geom_handle, &
      &                 key = 'geometry')
    ! Try to open a sub table of the geometry.
    call aot_table_open(L = conf, parent = geom_handle,  &
      &                 thandle = geom_subhandle, pos = 1)

    if (geom_handle .gt. 0) then
      ! Geometry table is defined
      if (geom_subhandle .eq. 0) then
        ! It is propably just a single entry, as the first entry is not a table
        ! in itself.
        call aot_table_close(L = conf, thandle = geom_subhandle) 
        call sdr_load_geom_table_single(me = me, conf = conf,  &
          &                             thandle = geom_handle, &
          &                             attr_pos = attr_pos,   &
          &                             transform = transform)
  
      else
        ! Geometry table contains propably multiple geometry definitions.
        call aot_table_close(L = conf, thandle = geom_subhandle)
        nTables = aot_table_length(L = conf, thandle = geom_handle)
        do iTab =1,nTables
          call aot_table_open(L = conf, parent = geom_handle, &
            &                 thandle = geom_subhandle, pos = iTab)

          call sdr_load_geom_table_single(me = me, conf = conf,     &
            &                             thandle = geom_subhandle, &
            &                             attr_pos = attr_pos,      &
            &                             transform = transform)

          call aot_table_close(L = conf, thandle = geom_subhandle)
        end do
      end if
    else
      write(logunit(0),*) 'ERROR loading geometry for spatial object: '
      write(logunit(0),*) '      NO geometry table defined!'
      call tem_abort()
    end if
   
    call aot_table_close(L = conf, thandle = geom_handle) 

  end subroutine sdr_load_geom_table
  ! *****************************************************************************
  

  ! *****************************************************************************
  !> This routine reads a single geometry table from the config file i.e 
  !! loading different geometry kinds like canoND, cube, periodic, STL etc.
  subroutine sdr_load_geom_table_single(me, conf, thandle, attr_pos, transform)
    ! --------------------------------------------------------------------------!
    type( sdr_geometry_type ), intent(inout) :: me
    type(flu_State) :: conf !< lua state
    integer, intent(in) :: thandle !< handle for spatialObject
    !> Position of the attribute to connect all primitives to.
    integer, intent(in) :: attr_pos
    type(tem_transformation_type), intent(in) :: transform
    ! --------------------------------------------------------------------------!
    integer :: iError
    character(len=labelLen) :: geom_kind
    ! --------------------------------------------------------------------------!

    ! Get the kind of the geometry to read.
    call aot_get_val(L = conf, thandle = thandle, key = 'kind', &
      &              val = geom_kind, ErrCode = iError)
    if (btest(iError, aoterr_Fatal)) then
      write(logunit(0),*) 'FATAL Error occured, while retrieving geometry kind:'
      if (btest(iError, aoterr_NonExistent)) &
        &  write(logunit(0),*) 'Geometry kind not defined!'
      if (btest(iError, aoterr_WrongType)) &
        & write(logunit(0),*) 'Geometry kind has wrong type!'
      call tem_abort()
    end if

    select case(upper_to_lower(trim(geom_kind)))
    case('canond', 'canonicalnd')
      call sdr_load_canonicalND(pntArray    = me%point,        &
        &                       lineArray   = me%line,         &
        &                       triArray    = me%triangle,     &
        &                       boxArray    = me%box,          &
        &                       spaObjArray = me%spatialObj,   &
        &                       attr_pos    = attr_pos,        &
        &                       transform   = transform,       &
        &                       conf = conf, thandle = thandle )
    case('triangle')
      call sdr_load_triangle(triArray    = me%triangle,   &
        &                    spaObjArray = me%spatialObj, &
        &                    attr_pos    = attr_pos,      &
        &                    transform   = transform,     &
        &                    conf        = conf,          &
        &                    thandle     = thandle        )
    case('stl')
      call sdr_load_stl(triArray    = me%triangle,   &
        &               spaObjArray = me%spatialObj, &
        &               attr_pos    = attr_pos,      &
        &               transform   = transform,     &
        &               conf        = conf,          &
        &               thandle     = thandle        )
    case('sphere')
      call sdr_load_sphere(sphArray    = me%sphere,     &
        &                  spaObjArray = me%spatialObj, &
        &                  attr_pos    = attr_pos,      &
        &                  transform   = transform,     &
        &                  conf        = conf,          &
        &                  thandle     = thandle        )
    case('ellipsoid')
      call sdr_load_ellipsoid(sphArray = me%ellipsoid,  &
        &                  spaObjArray = me%spatialObj, &
        &                  attr_pos    = attr_pos,      &
        &                  transform   = transform,     &
        &                  conf        = conf,          &
        &                  thandle     = thandle        )
    case('cylinder')
      call sdr_load_cylinder(cylArray  = me%cylinder,   &
        &                  spaObjArray = me%spatialObj, &
        &                  attr_pos    = attr_pos,      &
        &                  transform   = transform,     &
        &                  conf        = conf,          &
        &                  thandle     = thandle        )
    case('spacer')
      call sdr_load_spacer(sprInterwovenArray = me%spacerInterwoven, &
        &                  cylArray           = me%cylinder,         &
        &                  spaObjArray        = me%spatialObj,       &
        &                  attr_pos           = attr_pos,            &
        &                  transform          = transform,           &
        &                  conf               = conf,                &
        &                  thandle            = thandle              )
    case('periodic')
      call sdr_load_periodic(periArray   = me%periPlane,  &
        &                    spaObjArray = me%spatialObj, &
        &                    attr_pos    = attr_pos,      &
        &                    transform   = transform,     &
        &                    conf        = conf,          &
        &                    thandle     = thandle        )

    case default
      write(logunit(0),*) 'ERROR: Unknown geometry kind: '//trim(geom_kind)
      call tem_abort()
    end select

  end subroutine sdr_load_geom_table_single
  ! *****************************************************************************


  ! ****************************************************************************
  !> This routine created sphere objects and new attribute and extend 
  !! a list of spatial objects if node intersected boundary has distance 
  !! refine.
  subroutine sdr_append_distanceRefineObject(coord, dx, iLevel, geometry,      &
    &                   intersected_first, intersected_last, intersected_object)
    ! --------------------------------------------------------------------------!
    !> Coordinate of current node
    integer, intent(in) :: coord(4)
    real(kind=rk), intent(in) :: dx !< size of the element
    !> Current node level
    integer, intent(in) :: iLevel
    !> type which contains all geometry object infos
    type(sdr_geometry_type), intent(inout) :: geometry
    !> This node's first position in intersected_object
    integer, intent(in) :: intersected_first
    !> This node's last position in intersected_object
    integer, intent(in) :: intersected_last
    !> Growing array of intersected_objects
    type(grw_intArray_type), intent(in) :: intersected_object
    ! --------------------------------------------------------------------------!
    type(tem_sphere_type) :: loc_sphere
    type(sdr_spatialObj_type) :: spatialObj 
    type(sdr_attribute_type) :: attribute
    integer :: iDistRefine, iObject
    integer :: distRef_first, distRef_last
    integer :: obj_pos, attr_pos, uni_pos
    logical :: wasAdded
    logical :: attrChecked(geometry%attribute%dynArray%nVals)
    ! --------------------------------------------------------------------------!
    attrChecked = .false.

    ! Set origin of sphere, barcenter of current node
    loc_sphere%origin = geometry%universe%origin &
      &               + (real(coord(1:3), kind=rk) + 0.5_rk) * dx
    
    loc_sphere%only_surface = .false.

    ! Set attribute kind to refinement and label to distance
    attribute%kind = sdr_Refinement_object
    attribute%label = 'distance_refine'
    
    do iObject = intersected_first, intersected_last
      obj_pos = intersected_object%val(iObject)

      attr_pos = geometry%spatialObj%val(obj_pos)%attribute_position
      ! if attribute contains distance function and this attribute is
      ! not checked before
      if (geometry%attribute%dynArray%val(attr_pos)%distRefine_id > 0 &
        .and. .not. attrChecked(attr_pos) ) then
        attrChecked(attr_pos) = .true.
        distRef_first = geometry%attribute%dynArray%val(attr_pos)%distance_first
        distRef_last = geometry%attribute%dynArray%val(attr_pos)%distance_last
        do iDistRefine = distRef_first, distRef_last
          ! set sphere radius
          loc_sphere%radius = geometry%attribute%distRefine        &
            &                               %val(iDistRefine)%radius
          ! set attribute level  
          attribute%level = geometry%attribute%distRefine                 &
            &                                 %val(iDistRefine)%reach_level
          
          ! Create distance refine object only if reach_level>0
          ! and same as current level to reduce number of spatial objects
          if (attribute%level > 0 .and. attribute%level==iLevel) then
            ! append attribute to dynamic array of attribute
            call append(geometry%attribute%dynArray, attribute, &
              &         pos = attr_pos, &
              &         wasAdded = wasAdded)
  
            if (wasAdded) then
             ! Update the counters for refinement kind of object.
              call append(geometry%attribute%kindPos(attribute%kind), attr_pos)
              geometry%attribute%dynArray%val(attr_pos)%id &
                & = geometry%attribute%kindPos(attribute%kind)%nVals
  
              call append(geometry%attribute%uni_name(attribute%kind), &
                &         val = attribute%label, &
                          pos = uni_pos)
  
              geometry%attribute%dynArray%val(attr_pos)%uni_id = uni_pos
            end if     
  
            ! append sphere to sphereObject
            call append(geometry%sphere, loc_sphere)
  
            ! append spatial object
            spatialObj%attribute_position = attr_pos
            spatialObj%geometry_primitive = sphere
            spatialObj%primitive_position = geometry%sphere%nVals
  
            call append(geometry%spatialObj, spatialObj)
          end if  
        end do
      end if
    end do

  end subroutine sdr_append_distanceRefineObject
  ! ****************************************************************************


end module sdr_geometry_module
