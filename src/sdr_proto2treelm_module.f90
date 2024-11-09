! Copyright (c) 2012-2017, 2023 Harald Klimach <harald.klimach@dlr.de>
! Copyright (c) 2012-2015, 2019 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2012, 2014 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2012 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2015 Samuel Ziegenberg
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2016, 2018 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2018 Daniel Fleischer <daniel.fleischer@student.uni-siegen.de>
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
!> author: Kannan Masilamani
!! author: Jiaxing
!! Module to create the actual mesh out of the flooded protoTree.
!!
!! This module provides the functionality to refine the computational domain
!! down to the requested level everywhere and identify the boundary conditions
!! in all 26 directions of the elements, which have boundaries.
!!
module sdr_proto2treelm_module
  use env_module,            only: long_k, rk, newunit, isLittleEndian, &
    &                              globalMaxLevels, labelLen
  use tem_param_module,      only: qQQQ
  use tem_grow_array_module, only: grw_longArray_type, grw_intArray_type,    &
    &                              grw_realArray_type, grw_char2dArray_type, &
    &                              init, append, truncate, destroy
  use tem_dyn_array_module,  only: positionOfVal, append, dyn_intArray_type, &
    &                              init, destroy
  use tem_bc_prop_module,    only: tem_bc_prop_type
  use tem_propHead_module,   only: tem_propHead_type
  use tem_topology_module,   only: tem_directChildren, tem_firstIDatLevel, &
    &                              tem_CoordOfId, tem_IdOfCoord,           &
    &                              tem_parentOf, tem_LevelOf
  use tem_geometry_module,   only: tem_eligibleChildren, tem_BaryOfCoord,    &
    &                              tem_CoordOfReal, tem_BaryOfID, tem_ElemSize
  use tem_aux_module,        only: tem_abort
  use tem_logging_module,    only: logunit
  use tem_property_module,   only: prp_hasBnd, prp_solid, prp_hasQVal, &
    &                              prp_isColored, prp_hasPolynomial,   &
    &                              prp_noSolidification, prp_hasNormal
  use tem_global_module,     only: tem_global_type, dump_tem_global
  use treelmesh_module,      only: treelmesh_type


  use ply_oversample_module, only: ply_convertFromOversample
  use ply_poly_project_module, only: ply_poly_project_type,     &
    &                                ply_poly_project_fillbody, &
    &                                ply_poly_project_n2m
  use ply_prj_header_module, only: ply_prj_header_out
  use sdr_cube_module,       only: sdr_cube_type
  use sdr_prototree_module,  only: sdr_protoTree_type, levelValues_type, &
    &                              sdr_neighbor_in_proto
  use sdr_geometry_module,   only: sdr_geometry_type, is_intersecting
  use sdr_periodic_module,   only: sdr_periodicPlane_type
  use sdr_node_module,       only: isFlooded_bit, isLeaf_bit, isColored_bit, &
    &                              IntersectsBoundary_bit, isTarget_bit,     &
    &                              hasBoundary_bit, isFluidifyable_bit,      &
    &                              sdr_nodeColors, sdr_nodeProp_btest,       &
    &                              sdr_color_log2char, sdr_bitfieldColors,   &
    &                              isNoSolidification_bit
  use sdr_attribute_module,  only: sdr_Boundary_Object, sdr_Seed_Object
  use sdr_config_module,     only: sdr_confHead_type
  use sdr_boundary_module,   only: sdr_identify_boundary,    &
    &                              sdr_qVal_no_intersection, &
    &                              sdr_find_periodic_neighbor

  use aot_out_module,        only: aot_out_type, aot_out_val,            &
    &                              aot_out_open, aot_out_close,          &
    &                              aot_out_open_table, aot_out_close_table

  use sdr_timer_module,      only: timer_handle_proto2treelm, &
    &                              timer_handle_dumping_mesh
  use tem_timer_module,      only: tem_startTimer,tem_stopTimer

  implicit none

  !> Data type contains final fluid info dumped into mesh file
  type sdr_temData_type
    !> list of fluid treeIDs
    type(grw_longArray_type) :: treeID
    !> list of property bits for fluid nodes
    type(grw_longArray_type) :: propertyBits
    !> list of boundary iDs for fluid nodes in 26 directions
    type(grw_longArray_type) :: bc_ID(qQQQ)
    !> list of color characters
    type(grw_char2dArray_type) :: colors
    !> list of color characters indicating their subresolution status.
    type(grw_char2dArray_type) :: subres
    !> Count number of subresolved elements for each color.
    integer, allocatable :: color_subres_count(:)
    !> File unit for the subresolution information of each color.
    integer, allocatable :: color_subres_unit(:)
    !> Measure of the volume covered by each color
    real(kind=rk), allocatable :: color_volume(:)
    !> maxlevel in the fluid domain
    integer :: maxLevel
    !> minlevel in the fluid domain
    integer :: minLevel
    !> link-wise distances from boundary for fluid nodes in 26 directions
    type(grw_realArray_type) :: qVal(qQQQ)
    !> Surface normals in boundary elements for which it is to be stored
    type(grw_realArray_type) :: bc_normal(3)
    !> treelm mesh type contains only bounding cube info
    !! needed to identify boundary and compute q-values
    type( treelmesh_type ) :: meshUniverse
    !> number of fluidifyable ( solids )
    integer :: nSolids = 0
    !> Number of fluids on each level
    integer :: nFluids(globalMaxLevels)
    !> Projection workspace.
    type(ply_poly_project_type) :: projection
  end type sdr_temData_type


contains


  ! ****************************************************************************!
  !> This subroutine creates the treelmesh from the flooded prototree.
  !!
  !! It will create the treelmesh out of all flooded elements and refine those
  !! leaves intersected by a refinement object.
  !! This routine will also identify the elements with boundary conditions
  !! and store the boundary information accordingly.
  !! Please keep in mind, that the flooding only considered the 6 face
  !! neighbors, while we create boundary informations for all 26 neighbors.
  !! @todo Neighbors that do not intersect boundary objects.
  !! This means, that it might happen, that any of the 20 further neighbors
  !! does not actually intersect a boundary object. We need to deal with
  !! this case properly.
  !! For example, a boundary condition might be selected based on the
  !! adjacent face neighbors for the direction in question.
  subroutine sdr_proto2treelm(proto, geometry, temData, header)
    ! --------------------------------------------------------------------------!
    !> preliminary tree
    type(sdr_protoTree_type), intent(inout) :: proto
    !> contains all geometrical objects
    type(sdr_geometry_type), intent(in) :: geometry
    !> temData contains final mesh created by this routine
    type(sdr_temData_type), intent(inout) :: temData
    !> Header data.
    type(sdr_confHead_type), intent(inout) :: header
    ! --------------------------------------------------------------------------!
    integer :: minlevel
    integer :: ilevel
    integer :: iDir
    integer :: iColor
    type(levelValues_type) :: leVal
    ! --------------------------------------------------------------------------!

    call tem_startTimer( timerHandle = timer_handle_proto2treelm )


    !effective bounding info
    temData%meshUniverse%global%effOrigin = &
      & huge(temData%meshUniverse%global%effOrigin)
    temData%meshUniverse%global%effLength = &
      & -huge(temData%meshUniverse%global%effLength)

    minlevel = header%minlevel

    write(logunit(1),*) 'Creating treelmesh ...'

    ! Initialize the growing arrays for the new mesh data, they will have to
    ! contain at least as many elements, as we found flooded leaf nodes so far.
    call init(temData%treeID, length=proto%nFloodedLeaves)
    call init(temData%propertyBits, length=proto%nFloodedLeaves)

    temData%maxLevel = 0
    temData%minLevel = huge(temData%minLevel)
    temData%nFluids = 0

    do iDir=1,qQQQ
      call init(temData%bc_ID(iDir))
      call init(temData%qVal(iDir))
    end do

    call init(temData%bc_normal(1))
    call init(temData%bc_normal(2))
    call init(temData%bc_normal(3))

    call init( me     = temData%colors,        &
      &        width  = proto%node%nColorChars )

    call init( me     = temData%subres,        &
      &        width  = proto%node%nColorChars )
    allocate( temData%color_subres_count(proto%node%nColors) )
    allocate( temData%color_subres_unit(proto%node%nColors) )
    allocate( temData%color_volume(proto%node%nColors) )
    temData%color_volume = 0.0_rk
    temData%color_subres_count = 0

    if (header%subresolution%nLevels > 0) then
      ! If we are going to deal with subresolution, we need to fill the
      ! projection data now.
      call ply_poly_project_fillbody( me         = temData%projection,  &
        &                             proj_init  = header%subresolution &
        &                                                %projection,   &
        &                             scheme_dim = 3                    )
    end if

    ! Set the length of the level cube length to the complete bounding cube
    ! length for level=0.
    ! Set some auxilary data describing the current level.
    leVal%dx = geometry%universe%extent
    leVal%ID_offset = tem_FirstIdAtLevel(0) ! first treeID on this level
    leVal%level = 0 ! level count

    ! Start with the root node, and recursively iterate down to the leaves.
    ! As we follow the child ordering from the space filling curve, this
    ! recursive procedure will result in a sorted list of elements as required
    ! by treelm.
    call traverse_tree( node_pos     = 1,                   &
      &                 leVal        = leVal,               &
      &                 proto        = proto,               &
      &                 geometry     = geometry,            &
      &                 temData      = temData,             &
      &                 header       = header,              &
      &                 meshUniverse = temData%meshUniverse )

    ! Close subelement color files again.
    do iColor=1,proto%node%nColors
      if (temData%color_subres_count(iColor) > 0) then
        close(temData%color_subres_unit(iColor))
      end if
    end do

    write(logunit(1),*) '                      done'
    write(logunit(2),*) 'Total number of fluids: ', temData%treeID%nVals
    write(logunit(2),*) ''

    if (temdata%minlevel < temdata%maxlevel) then
      ! For multilevel meshes print a summary on element distribution:
      write(logunit(2),'(2(a,i0))') 'Distributed across levels ', &
        &                           temdata%minlevel, ' to ', temdata%maxlevel
      do iLevel=temdata%minlevel,temdata%maxlevel
        write(logunit(2),'(2(a,i0),a)') 'Level ', iLevel, ': ', &
          &                             temdata%nFluids(iLevel), ' Elements'
      end do
      write(logunit(2),*) ''
    end if

  call tem_stopTimer( timerHandle = timer_handle_proto2treelm )


  end subroutine sdr_proto2treelm
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> Recursively traverse the tree in a depth first manner to obtain the
  !! the ordering required by treelm.
  !!
  !! Only consider flooded nodes here, no need to traverse further down any
  !! non-flooded node, as it does not belong to the computational domain.
  !! All virtual nodes containing at least one flooded leaf node has already
  !! been determined and set after flooding in sdr_flood_module::sdr_flood_tree.
  !! If the leaf node is reached, hand over to the ::refine_leaf routine, to
  !! find all elements, that actually should be created in the computational
  !! domain.
  !!
  !! @todo KM: Store has_boundary information for every leaf node in
  !! sdr_refine_leaf routine to avoid excessive boundary checking far
  !! away from the boundaries.
  recursive subroutine traverse_tree(node_pos, leVal, proto, geometry, &
    &                                temData, header, meshUniverse)
    ! ---------------------------------------------------------------------- !
    !> Position of leaf in the preliminary tree
    integer, intent(in) :: node_pos
    !> level value of current node
    type(levelValues_type), intent(in) :: leVal
    !> preliminary tree
    type(sdr_protoTree_type), intent(in) :: proto
    !> contains all geometrical objects
    type(sdr_geometry_type), intent(in) :: geometry
    !> temData contains final treelmesh
    type(sdr_temData_type), intent(inout) :: temData
    !> Header information
    type(sdr_confHead_type), intent(in) :: header
    !> contains bounding cube information
    type(treelmesh_type), intent(inout) :: meshUniverse
    ! ---------------------------------------------------------------------- !
    integer :: minlevel
    integer :: child_pos
    integer :: iChild
    logical :: check_bnd
    type(levelValues_type) :: child_leVal
    integer :: nodeprops
    character :: nodecolors(proto%node%nColorChars)
    ! ---------------------------------------------------------------------- !

    minlevel = header%minlevel

    nodeprops = ibits(proto%node%PropertyBits &
      &                         %val(proto%node%propLength, Node_pos), &
      &               pos = proto%node%lastbyte_pos, len=8)

    !> If this node is flooded or fluidifyable it needs to be included in the
    !! final mesh.
    if ( btest(nodeprops, isFlooded_bit) ) then

      ! If flooded take some action, depending on the node type:
      ifleaf: if ( btest(nodeprops, isLeaf_bit) ) then
        ! This node IS A LEAF!

        ! Get its colors:
        nodecolors = sdr_nodeColors(proto%node, Node_pos)

        ! If this node is colored, set its bit accordingly.
        if (any(ichar(nodecolors) > 0)) then
          nodeprops = ibset(nodeprops, isColored_bit)
        end if

        ! Add this leaf or if necessary its children down to the configured
        ! resolution to the list of fluid elements and find the correct boundary
        ! conditions.
        ! JQ: no need to check boundary condition for isFluidifyable nodes
        ! Check boundary only if it is marked with hasBoundary_bit
        check_bnd = ( .not. btest(nodeprops, isFluidifyable_bit) .and. &
          &           btest(nodeprops, hasBoundary_bit) )

        ! No need for further refinement.
        ! Just append the node to temData
        call proto2Treelm(                                       &
          &    node_pos       = node_pos,                        &
          &    treeID         = proto%node%treeID%val(node_pos), &
          &    nodeprops      = nodeprops,                       &
          &    nodecolors     = nodecolors,                      &
          &    leVal          = leVal,                           &
          &    proto          = proto,                           &
          &    geometry       = geometry,                        &
          &    check_bnd      = check_bnd,                       &
          &    temData        = temData,                         &
          &    meshUniverse   = meshUniverse                     )

      else ifleaf

        ! NOT A LEAF
        ! Check if it is a target element, if so add it as an element and create
        ! subelement resolution data.
        ! Otherwise traverse further down through all eight children.
        iftarget: if ( btest(nodeprops, isTarget_bit) ) then
          ! Target nodes are to be added to the final treelmesh.
          ! The children below them are only used to describe subelement
          ! resolved boundaries with the help of polynomials.

          call create_target( node_pos     = node_pos,            &
            &                 proto        = proto,               &
            &                 geometry     = geometry,            &
            &                 leVal        = leVal,               &
            &                 meshUniverse = meshUniverse,        &
            &                 header       = header,              &
            &                 temData      = temData              )

        else iftarget

          ! This is a virtual element.
          ! Set some auxilary data describing the next level (common to all
          ! children).
          child_leVal%level = leVal%level + 1
          child_leVal%dx = 0.5_rk * leVal%dx
          child_leVal%ID_offset = tem_FirstIdAtLevel( child_leVal%level )
          do iChild = 0, 7
            child_pos = proto%node%linkpos(1)%val(node_pos) + iChild
            call traverse_tree( node_pos     = child_pos,   &
              &                 leVal        = child_leVal, &
              &                 proto        = proto,       &
              &                 geometry     = geometry,    &
              &                 temData      = temData,     &
              &                 header       = header,      &
              &                 meshUniverse = meshUniverse )
          end do

        end if iftarget

      end if ifleaf
    end if !if flooded

  end subroutine traverse_tree
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> Routine to convert protoTree to Treelm data format.
  !! append all leaves to the temData%treeID
  !!
  !! Add this element to the list of elements in the final tree, by adding:
  !!       * treeID
  !!       * propertyBits
  !! If there is a boundary, also add the correct boundary data into the
  !! 26 direct neighbors (bc_ID), append to all 26 even if there is no
  !! boundary condition in the given direction, no boundaries are indicated
  !! by 0, boundary conditions are decided based on the trumping rule, that
  !! is the minimal bc_ID of all boundary objects in the corresponding
  !! neighbor are chosen.
  !! Treat periodic boundaries here if possible by putting the opposite
  !! treeID into the bc_ID.
  subroutine proto2Treelm( node_pos, treeID, leVal, proto, geometry, check_bnd,&
    &                      temData, nodeprops, nodecolors, meshUniverse )
    ! ---------------------------------------------------------------------- !
    !> Position of leaf in the preliminary tree
    integer, intent(in) :: node_pos
    !> treeID of current node
    integer(kind=long_k), intent(in) :: treeID
    !> propertyBits of current node
    integer, intent(in) :: nodeprops
    !> Color information of the cuurent node
    character, intent(in) :: nodecolors(:)
    !> level value of current node
    type(levelValues_type), intent(in) :: leVal
    !> preliminary tree
    type(sdr_protoTree_type), intent(in) :: proto
    !> contains all geometrical objects
    type(sdr_geometry_type), intent(in) :: geometry
    !> does this node has boundary neighbor
    logical, intent(in) :: check_bnd
    type(sdr_temData_type), intent(inout) :: temData
    !> contains bounding cube information
    type(treelmesh_type), intent(inout) :: meshUniverse
    ! ---------------------------------------------------------------------- !
    integer :: iDir
    integer(kind=long_k) :: BC_ID(qQQQ)
    real(kind=rk) :: qVal( qQQQ ) !< link-wise distance for boundary elements
    real(kind=rk) :: bc_normal(3) ! normal of closest surface
    integer(kind=long_k) :: propBits
    type(sdr_cube_type) :: node_cube
    integer :: node_coord(4)
    real(kind=rk) :: elemvol
    integer :: nodepropbits(proto%node%propLength)
    integer :: iColor, col_int, col_bit
    ! ---------------------------------------------------------------------- !

    ! Intitalize boundary info.
    BC_ID = 0
    qVal = sdr_qVal_no_intersection
    ! coordinate of current node treeID
    node_coord = tem_coordOfId(treeID = treeID,         &
      &                        offset = leVal%ID_offset )
    ! Find which neighbors have boundaries.
    ! save bcID in `BC_ID`, q-value in `qVal` and surface normal in `bc_normal`
    if (check_bnd) then
      call sdr_identify_boundary(node_pos, treeID, node_coord, leVal, proto, &
        &                        geometry, BC_ID, qVal, bc_normal,           &
        &                        meshUniverse                                )
    end if

    ! Define the cube to intersect with
    node_cube%origin = leVal%dx * node_coord(1:3) &
      &              + geometry%universe%origin
    node_cube%halfwidth = 0.5_rk*leVal%dx
    node_cube%center = node_cube%origin + node_cube%halfwidth
    node_cube%extent = leVal%dx

    elemvol = leval%dx**3

    ! update maxLevel and minLevel according to current level
    temdata%maxlevel = max(temdata%maxlevel, leVal%level)
    temdata%minlevel = min(temdata%minlevel, leVal%level)
    temdata%nFluids(leval%level) = temdata%nFluids(leval%level) + 1

    ! Append this node to the temData lists
    call append( temData%treeID, treeID )

    ! compute effective bounding cube origin and length using node origin
    ! and its extent+origin respectively
    meshUniverse%global%effOrigin &
      & = min(meshUniverse%global%effOrigin, node_cube%origin)
    meshUniverse%global%effLength &
      & = max(meshUniverse%global%effLength, node_cube%origin + leVal%dx)

    propBits = 0_long_k

    ! The color of the current node.
    if ( btest(nodeprops, isColored_bit) ) then
      propBits = ibset(propBits, prp_isColored)
      call append(temData%colors, nodecolors)
      ! Check the colored status, to account for the volume fraction of each
      ! color.
      nodepropbits = proto%node%PropertyBits%val(:,Node_pos)
      do iColor=1,proto%node%nColors
        ! Figure out the bit 0 in the integer field for this color
        col_int = (iColor-1) / proto%node%bytes_per_int + 1
        col_bit = mod(iColor-1, proto%node%bytes_per_int)*8
        if ( btest(nodepropbits(col_int), col_bit+1) ) then
          temData%color_volume(iColor) = temData%color_volume(iColor) &
            &                          + elemvol
        end if
      end do
    end if

    ! Decide on hasBnd, and set temData%PropertyBits and temData%BC_ID
    ! accordingly.
    if ( any(BC_ID /= 0_long_k) ) then
      ! there is an adjacent boundary
      propBits = ibset( propBits, prp_hasBnd )

      ! JQ: if it has qVal, set propertyBits and append qVal
      if ( any( qVal > 0.0_rk )) then
        propBits = ibset( propBits, prp_hasQVal )
        do iDir = 1,qQQQ
          call append( temData%qVal(iDir), qVal(iDir) )
        end do
      end if

      if ((bc_normal(1)**2 + bc_normal(2)**2 + bc_normal(3)**2) > 0.0_rk) then
        propBits = ibset(propBits, prp_hasNormal)
        do iDir = 1,3
          call append(temData%bc_normal(iDir), bc_normal(iDir))
        end do
      end if

      ! append bc_id
      do iDir = 1,qQQQ
        call append( temData%BC_ID(iDir), BC_ID(iDir) )
      end do

    end if

    if ( btest(nodeprops, isFluidifyable_bit) ) then
      propBits = ibset(propBits, prp_solid)
      temData%nSolids = temData%nSolids + 1
    end if

    ! If this node is a noSolidification node, set its bit accordingly.
    if ( btest(nodeprops, isNoSolidification_bit) ) then
      propbits = ibset(propbits, prp_noSolidification)
    end if

    ! append propertyBits
    call append( temData%PropertyBits, propBits )

  end subroutine proto2Treelm
  ! ****************************************************************************!


  ! ****************************************************************************!
  !> Routine to create an element with subelement resolution.
  !!
  !! This routine is used to define a target element, that is a node, that
  !! reached its final level, but contains boundaries, that are to be resolved
  !! on a subelement basis. For the treelmesh, this is just a normal element
  !! with additional information attached to it. That additional information
  !! is covered by a property.
  !! For the actual information in the subelement resolution, the children of
  !! the target node need to be considered.
  !! We use a separate routine for this case, as it is ensured, that this
  !! element never will be refined, due to adjacent boundaries. (Any such
  !! refinement was already done during the building of the protoTree).
  !! Also we will ignore Q-Value complications here for the identification of
  !! boundaries. If subelement resolution is active, no Q-Values might be used
  !! and the other way around.
  subroutine create_target( node_pos, proto, geometry, leVal, &
    &                       meshUniverse, header, temData )
    ! --------------------------------------------------------------------------!
    !> Position of leaf in the preliminary tree
    integer, intent(in) :: node_pos
    !> preliminary tree
    type(sdr_protoTree_type), intent(in) :: proto
    !> contains all geometrical objects
    type(sdr_geometry_type), intent(in) :: geometry
    !> level value of parent node
    type(levelValues_type), intent(in) :: leVal
    !> contains bounding cube information
    type(treelmesh_type), intent(inout) :: meshUniverse
    !> Configuration header.
    type(sdr_confHead_type), intent(in) :: header
    !> Data for the final mesh.
    type(sdr_temData_type), intent(inout) :: temData
    ! --------------------------------------------------------------------------!
    ! Colors, that are found within this element.
    character :: nodecolors(proto%node%nColorChars)
    character :: subresolved_colors(proto%node%nColorChars)
    integer :: nodepropbits(proto%node%propLength)

    logical :: color_subresolved(proto%node%nColors)
    logical :: color_subresflooded(proto%node%nColors)
    integer(kind=long_k) :: nodeID, propBits
    integer :: nObjects
    integer :: obj_pos
    integer :: iDir, iObject, iColor
    integer :: attr_kind, attr_pos
    integer :: bc_color, bcid
    integer :: col_int, col_bit
    integer :: node_coord(4)
    type(sdr_cube_type) :: node_cube
    integer(kind=long_k) :: BC_ID(qQQQ)
    integer(kind=long_k) :: fallback_bc
    integer :: neighbor_pos, neighbor_level
    real(kind=rk) :: dummy_qval
    logical :: testall
    character(len=LabelLen) :: col_name
    character(len=4) :: EndianSuffix
    integer :: subres_pos
    real(kind=rk), allocatable :: color_modes(:,:)
    real(kind=rk), allocatable :: color_oversampled(:,:)
    real(kind=rk), allocatable :: color_nodes(:,:)
    real(kind=rk) :: elemvol
    integer :: rl
    integer :: nDofs
    integer :: nPoints
    integer :: nOversampled
    integer :: record
    integer :: intersected_first, intersected_last
    logical :: unKnownBnd
    ! --------------------------------------------------------------------------!

    if (isLittleEndian) then
      EndianSuffix = '.lsb'
    else
      EndianSuffix = '.msb'
    end if

    propBits = 0_long_k

    nodeID = proto%node%treeID%val(node_pos)

    ! If no boundary is found in a direction, facing a non-flooded element, we
    ! fallback to the minimal boundary condition in the target element itself.
    fallback_bc = int(proto%node%minBCID%val(node_pos), kind=long_k)

    if (fallback_bc == -1_long_k) then
      ! Periodic boundaries not allowed in subresolved elements !
      write(logunit(1),'(a)') 'There is a periodic boundary running through a'
      write(logunit(1),'(a,i0,a)') 'subresolved element (treeID=', nodeID, &
        &                          ').'
      write(logunit(1),'(a)') 'This can not be properly handled!'
      write(logunit(1),'(a)') 'Please ensure that periodic boundaries and'
      write(logunit(1),'(a)') 'subelement resolved boundaries are sufficiently'
      write(logunit(1),'(a)') 'far apart!'
      write(logunit(1),*)
      write(logunit(1),'(a)') 'STOPPING!'
      call tem_abort()
    end if

    ! update maxLevel and minLevel according to current level
    temdata%maxlevel = max(temdata%maxlevel, leVal%level)
    temdata%minlevel = min(temdata%minlevel, leVal%level)

    ! Define the intersected cube
    node_coord = tem_coordOfId( treeID = nodeID,         &
      &                         offset = leVal%ID_offset )
    node_cube%origin = leVal%dx * node_coord(1:3) &
      &              + geometry%universe%origin
    node_cube%halfwidth = 0.5_rk*leVal%dx
    node_cube%center = node_cube%origin + node_cube%halfwidth
    node_cube%extent = leVal%dx

    elemvol = leVal%dx**3

    ! compute effective bounding cube origin and length using node origin
    ! and its extent+origin respectively
    meshUniverse%global%effOrigin &
      &  = min(meshUniverse%global%effOrigin, node_cube%origin)
    meshUniverse%global%effLength &
      &  = max(meshUniverse%global%effLength, node_cube%origin + leVal%dx)

    ! Look up neighbors and set boundaries accordingly
    do iDir=1,qQQQ

      BC_ID(iDir) = 0_long_k ! Default to no boundary.
      unKnownBnd = .false.
      neighbor_pos = sdr_neighbor_in_proto( proto, node_coord, iDir,  &
        &                                   neighbor_level )

      ! Only if the neighbor was not flooded, we need a boundary.
      ! Flooded neighbors might be target elements or plain fluids, if there
      ! is a boundary, but no target, the element will not be flooded.
      if ( .not. sdr_nodeProp_btest(node  = proto%node,   &
        &                           iNode = neighbor_pos, &
        &                           bit   = isFlooded_bit ) ) then
        ! The node is not flooded, if it is intersected, we can look up the
        ! boundary ID to use in this direction.
        if (sdr_nodeProp_btest(node  = proto%node,           &
          &                    iNode = neighbor_pos,         &
          &                    bit   = intersectsBoundary_bit) ) then
          BC_ID(iDir) = int(proto%node%minBCID%val(neighbor_pos), kind=long_k)
          ! Treating periodic boundaries identified by an ID of -1:
          if (BC_ID(iDir) == -1_long_k) then
            call sdr_find_periodic_neighbor( elemBary       = node_cube%center,&
              &                              iDir           = iDir,            &
              &                              bc_id          = BC_ID(iDir),     &
              &                              qVal           = dummy_qVal,      &
              &                              unKnownBnd     = unKnownBnd,      &
              &                              neighbor_pos   = neighbor_pos,    &
              &                              neighbor_level = neighbor_level,  &
              &                              leVal          = leVal,           &
              &                              proto          = proto,           &
              &                              geometry       = geometry,        &
              &                              meshUniverse   = meshUniverse     )
          end if
        else
          ! No boundary in this direction, but no flooded element, use the
          ! fallback boundary.
          BC_ID(iDir) = fallback_bc
        end if
        ! if periodic boundary encounters hanging node with no property then
        ! use fall back boundary
        if (unKnownBnd) BC_ID(iDir) = fallback_bc
      end if

    end do

    if ( any(BC_ID /= 0_long_k) ) then
      ! there is at least one adjacent boundary
      propBits = ibset( propBits, prp_hasBnd )
      ! append bc_id
      do iDir = 1,qQQQ
        call append( temData%BC_ID(iDir), BC_ID(iDir) )
      end do
    end if

    ! test all spatial objects only for root node.
    ! Check whether to intersect with all geometries
    testAll = (leVal%level == 0)

    if (testAll) then
      nObjects = geometry%spatialObj%nVals
      intersected_first = 1
      intersected_last = nObjects
    else
      intersected_first = proto%node%userObjPos%val(node_pos)%first
      intersected_last = proto%node%userObjPos%val(node_pos)%last
      nObjects = intersected_last - intersected_first + 1
    end if

    ! Check which colors actually need to be subresolved.
    color_subresolved = .false.

    objLoop: do iObject = intersected_first, intersected_last

      if (testAll) then
        obj_pos = iObject
      else
        obj_pos = proto%node%intersected_object%val(iObject)
      end if

      attr_pos = geometry%spatialObj%val(obj_pos)%attribute_position
      attr_kind = geometry%attribute%dynArray%val(attr_pos)%kind

      if (attr_kind == sdr_Boundary_object) then
        bcid = geometry%attribute%dynArray%val(attr_pos)%id
        bc_color = geometry%attribute%bc_color_id(bcid)
        if (bc_color > 0) then
          color_subresolved(bc_color) = color_subresolved(bc_color)       &
            &                           .or. geometry%attribute%dynArray  &
            &                                %val(attr_pos)%subresolution
        end if
      end if

    end do objLoop

    ! In retrospect, we now have to remove the flooding status for those colors
    ! with an intersecting boundary but no subresolution.
    nodepropbits = proto%node%PropertyBits%val(:,Node_pos)
    color_subresflooded = .false.
    do iColor=1,proto%node%nColors
      ! Figure out the bit 0 in the integer field for this color
      col_int = (iColor-1) / proto%node%bytes_per_int + 1
      col_bit = mod(iColor-1, proto%node%bytes_per_int)*8

      if ( btest(nodepropbits(col_int), col_bit) ) then
        ! If this color is intersected, clear the flooded status, if it is not
        ! also to be subresolved.
        if ( .not. color_subresolved(iColor) ) then
          ! Bit 1 for this color indicates the flooding status, clear it if
          ! there is an intersection, but no subresolution.
          nodepropbits = ibclr(nodepropbits(col_int), col_bit+1)
        end if
      end if

      ! If the element is flooded by the color, we set the subresolution action
      ! depending on the intersection status. Ignore the none color here.
      if (iColor /= proto%node%none_color_id) then
        if ( btest(nodepropbits(col_int), col_bit+1) ) then
          ! If the element is flooded by any color, set the isColored bit.
          propBits = ibset(propBits, prp_isColored)

          ! Check each color for its intersection status
          ! This will only be true, if the element is flooded by the color and
          ! intersected by a boundary of this color with the subresolution flag.
          color_subresflooded(iColor) = btest(nodepropbits(col_int), col_bit)
          ! Count fully flooded elements for the volume fraction of the color.
          if (.not. color_subresflooded(iColor)) then
            temData%color_volume(iColor) = temData%color_volume(iColor) &
              &                          + elemvol
          end if
        end if
      end if
    end do

    if ( btest(propBits, prp_isColored)) then
      nodecolors = sdr_bitfieldColors( node     = proto%node,  &
        &                              bitfield = nodepropbits )
      call append(temData%colors, nodecolors)
    end if

    if (any(color_subresflooded)) then
      ! This element provides polynomial informations on the colors, thus
      ! hasPolynomial bit also needs to be set.
      propBits = ibset(propBits, prp_hasPolynomial)
      ! Append the bitmask on, which colors are to be found in this element to
      ! the list of colors.
      subresolved_colors = sdr_color_log2char( node     = proto%node,       &
        &                                      logicals = color_subresolved )
      call append(temData%subres, subresolved_colors)
    end if

    ! No need for further refinement.
    ! Append this node to the temData lists
    call append( temData%treeID, nodeID )
    call append( temData%PropertyBits, propBits )


    nDofs = temData%projection%body_3d%ndofs
    noversampled = (temData%projection%oversamp_degree+1)**3
    nPoints = temData%projection%nQuadPointsPerDir**3

    allocate(color_modes(nDofs,1))
    allocate(color_oversampled(noversampled,1))
    allocate(color_nodes(nPoints,1))

    inquire(iolength = rl) color_modes

    do iColor=1,proto%node%nColors
      ! Iterate over all colors that are to be subresolved and actually flood
      ! this element. Only for these we need to create polynomial data.
      if ( color_subresflooded(iColor) ) then
        col_name = geometry%attribute%uni_name(sdr_Seed_object)%val(iColor)
        subres_pos = positionofval( header%subresolution%color_values%label, &
          &                         trim(col_name)                           )
        write(logunit(6),*) 'Creating subresolution data for color ' &
          &                 //trim(col_name)
        write(logunit(6),*) 'in element ', nodeID

        temData%color_subres_count(iColor) = temData &
          &                                  %color_subres_count(iColor) + 1
        record = temData%color_subres_count(iColor)

        if (temData%color_subres_count(iColor) == 1) then
          ! First subresolution info for this color, create a file for it
          temData%color_subres_unit(iColor) = newunit()
          open( unit = temData%color_subres_unit(iColor),                  &
            &   file = trim(header%folder) // 'subresdata_'                &
            &          // trim(col_name) // EndianSuffix,                  &
            &   action = 'write', access = 'direct', form = 'unformatted', &
            &   recl = rl, status = 'replace'                              )
        end if

        ! Forall color_nodes: sdr_point_color
        ! nodal2modal(color_nodes, color_modes)
        ! Maybe do this inplace (without separated color_nodes)?
        call sdr_color_points(                                              &
          &       nPoints      = nPoints,                                   &
          &       nodals       = color_nodes(:,1),                          &
          &       point        = temData%projection%body_3d%nodes,          &
          &       fill         = header%subresolution%color_values          &
          &                            %fill%val(subres_pos),               &
          &       void         = header%subresolution%color_values          &
          &                            %void%val(subres_pos),               &
          &       target_pos   = node_pos,                                  &
          &       proto        = proto,                                     &
          &       iColor       = iColor                                     )

        call ply_poly_project_n2m( me         = temData%projection, &
          &                        dim        = 3,                  &
          &                        nVars      = 1,                  &
          &                        nodal_data = color_nodes,        &
          &                        modal_data = color_oversampled   )

        call ply_convertFromOversample(modalCoeffs = color_oversampled,  &
          &                            poly_proj   = temData%projection, &
          &                            nDim        = 3,                  &
          &                            state       = color_modes         )

        ! Account the volume fraction for the color.
        temData%color_volume(iColor) = temData%color_volume(iColor) &
          &                          + elemvol * color_modes(1,1)
        write(temData%color_subres_unit(iColor), rec=record) color_modes
      end if
    end do

    deallocate(color_nodes)
    deallocate(color_oversampled)
    deallocate(color_modes)

  end subroutine create_target
  ! ****************************************************************************!



  ! ****************************************************************************
  !> Get the color at all given points.
  !!
  !! Points need to be in the interval [-1,1].
  subroutine sdr_color_points( nodals, nPoints, point, target_pos, fill, void, &
    &                          proto, iColor )
    ! --------------------------------------------------------------------------!
    real(kind=rk), intent(out) :: nodals(:)
    integer, intent(in) :: nPoints
    !> Point for which to evaluate the color.
    real(kind=rk), intent(in) :: point(:,:)
    !> Position of leaf in the preliminary tree
    integer, intent(in) :: target_pos
    real(kind=rk), intent(in) :: fill, void
    !> preliminary tree
    type(sdr_protoTree_type), intent(in) :: proto
    integer, intent(in) :: iColor
    ! --------------------------------------------------------------------------!
    integer :: subelem
    integer :: iPoint
    integer :: col_int, col_bit
    integer :: iChild
    real(kind=rk) :: child_origin(3)
    real(kind=rk) :: child_point(3)
    integer :: point_coord(3)
    ! --------------------------------------------------------------------------!

    col_int = (iColor-1) / proto%node%bytes_per_int + 1
    col_bit = mod(iColor-1, proto%node%bytes_per_int)*8 + 1

    do iPoint=1,nPoints

      ! Switch to next point in list
      child_point = point(iPoint,:)

      ! Reset to target element
      subelem = target_pos

      ! Loop over levels (until value is found)
      lvl_loop: do

        ! Polynomial points in the range (-1,1), coord(point) = ceil(point)
        ! Use maximum to ensure integer coordinates in {0,1}
        point_coord = max(ceiling(child_point),0)
        iChild = point_coord(1) + 2*point_coord(2) + 4*point_coord(3)
        subelem = proto%node%linkpos(1)%val(subelem) + iChild

        ! Rescale point position to the child element
        child_origin = [-0.5_rk, -0.5_rk, -0.5_rk] + point_coord
        child_point = (child_point - child_origin)*2

        if ( btest(proto%node%PropertyBits%val(col_int, subelem), &
          &  col_bit) ) then
          ! Flooded, check wether we reached the leaf
          if ( sdr_nodeProp_btest( node  = proto%node, &
            &                      iNode = subelem,    &
            &                      bit   = isLeaf_bit  ) ) then
            ! Reached a flooded leaf, return the foreground
            nodals(iPoint) = fill
            EXIT lvl_loop
          end if
        else
          ! Not flooded: done return background
          nodals(iPoint) = void
          EXIT lvl_loop
        end if

      end do lvl_loop

    end do

  end subroutine sdr_color_points
  ! ****************************************************************************!



  ! ****************************************************************************
  !> This routine dumps the final fluid tree leaves in the disk
  subroutine sdr_dump_treelm( temData, geometry, confHeader )
    ! --------------------------------------------------------------------------!
    !> temData contains final treelmesh
    type(sdr_temData_type), intent(in) :: temData
    !> Bounding cube, the prototree lives in.
    type(sdr_geometry_type), intent(in) :: geometry
    !> config header info.
    type(sdr_confHead_type), intent(inout) :: confHeader
    ! --------------------------------------------------------------------------!
    character(len=4) :: EndianSuffix
    type(tem_global_type) :: tem_global
    integer :: iElem
    integer :: meshunit, bndunit, rl, distunit, colunit
    character(len=300) :: ElemFileName
    character(len=300) :: BndFileName
    character(len=300) :: ColFileName
    character(len=300) :: BndDistFileName  !< Boundary Distance File Name
    character(len=LabelLen) :: col_name
    type(aot_out_type) :: conf !< aotus lua state to write output
    integer(kind=long_k) :: BC_id(qQQQ)
    real(kind=rk) :: qVal(qQQQ)
    real(kind=rk) :: bc_normal(3)
    integer :: iDir, iAttr, iProp, iBC, iColor
    integer :: nProperties_max
    type(tem_propHead_type), allocatable :: property(:)
    integer :: nColors
    integer :: nBCs
    integer :: vpos
    integer(kind=long_k) :: id
    integer :: uni_id
    integer :: iUnibc
    real(kind=rk) :: cval
    type(dyn_intArray_type) :: actBC
    ! --------------------------------------------------------------------------!

    call tem_startTimer( timerHandle = timer_handle_dumping_mesh )

    write(logunit(1),*) 'dumping treelmesh to dir: '//trim(confHeader%folder)

    if (isLittleEndian) then
      EndianSuffix = '.lsb'
    else
      EndianSuffix = '.msb'
    end if

    ! Initialize nProperties
    tem_global%nProperties = 0
    nProperties_max = 4

    ! allocate property list and propHeader infos
    allocate( property( nProperties_max ) )

    ! ---------------------------------------------------------------------------
    ! Dump boundary only if bnd nElems > 0
    ! ---------------------------------------------------------------------------
    ! Dump boundary info to bnd.lua
    if (temData%bc_ID(1)%nVals > 0) then

      tem_global%nProperties = tem_global%nProperties + 1
      ! property 1 is boundary
      property(tem_global%nProperties)%bitpos = prp_hasBnd
      property(tem_global%nProperties)%label  = 'has boundaries'
      property(tem_global%nProperties)%nElems = temData%bc_ID(1)%nVals

      ! -------------------------------------------------------------------------
      ! Dump boundary ID to bnd.lsb or bnd.msb
      bndunit = newunit()
      BndFileName = trim(confHeader%folder) // 'bnd' // EndianSuffix

      inquire(iolength = rl) BC_id(:)
      open(unit = bndunit, file = trim(BndFileName), action = 'write', &
        &  form = 'unformatted', access = 'direct', &
        &  recl = rl)

      do iElem=1,temData%BC_ID(1)%nVals
        do iDir=1,qQQQ
          ! Boundaries might be defined, but not actually appear in the mesh.
          ! This is especially true for multi-color meshes, where there can be
          ! interfaces between colors, without a real color.
          ! We sort out these boundaries and renumber the boundary references
          ! accordingly. Negative values, that indicate other elements, need
          ! to be maintained.
          ! The mapping from the final bc_ids to the original unique boundaries
          ! is stored in actBC.
          id = temData%BC_ID(iDir)%val(iElem)
          if ( id > 0_long_k ) then
            iAttr = geometry%attribute%kindpos(sdr_Boundary_object) &
              &                       %val(id)
            uni_id = geometry%attribute%dynArray%val(iAttr)%uni_id
            call append( me  = actBC,  &
              &          val = uni_id, &
              &          pos = iBC     )
            BC_id(iDir) = int(iBC, kind=long_k)
          else
            BC_id(iDir) = temData%BC_ID(iDir)%val(iElem)
          end if
        end do
        write(bndunit, rec=iElem) BC_id
      end do

      nBCs = actBC%nVals
      call aot_out_open( conf, trim(confHeader%folder)//'bnd.lua' )
      call aot_out_val( put_conf = conf, vname = 'nSides', val = qQQQ )
      call aot_out_val( put_conf = conf, vname = 'nBCtypes', &
        &               val = nBCs )

      call aot_out_open_table( conf, 'bclabel' )
      do iBC = 1,nBCs
        iUnibc = actBC%val(iBC)
        call aot_out_val( conf,                                        &
          &               val = geometry%attribute                     &
          &                             %uni_name(sdr_Boundary_object) &
          &                             %val(iUnibc)                   )
      end do
      call aot_out_close_table(conf)
      call aot_out_close(conf)

      close(bndunit)
    end if
    ! ---------------------------------------------------------------------------

    ! ---------------------------------------------------------------------------
    ! Dump qVal only if qVal nElems > 0
    ! ---------------------------------------------------------------------------
    ! JQ: Dump qVal info to qval.lua
    if (temData%qVal(1)%nVals > 0) then
      tem_global%nProperties = tem_global%nProperties + 1
      ! JQ: add prp_hasQVal property information
      ! propery 2 is has qval
      property(tem_global%nProperties)%bitpos = prp_hasQVal
      property(tem_global%nProperties)%label  = 'has qVal'
      property(tem_global%nProperties)%nElems = temData%qVal(1)%nVals
      call aot_out_open( conf, trim(confHeader%folder)//'qval.lua' )

      ! Dump logical variables that indicate whether the corresponding BC
      ! has qVal
      call aot_out_open_table( conf, 'hasQVal' )
      do iBC = 1, actBC%nVals
        iUnibc = actBC%val(iBC)
        call aot_out_val( conf, &
          &               val = geometry%attribute%bc_uni_calcdist%val(iUnibc) )
      end do
      call aot_out_close_table(conf)
      call aot_out_close(conf)
      ! -------------------------------------------------------------------------

      ! -------------------------------------------------------------------------
      ! JQ: Dump boundary qVal to file qval.lsb
      distUnit = newunit()
      BndDistFileName = trim(confHeader%folder) // 'qval' // EndianSuffix

      inquire(iolength = rl) qVal(:)
      open(unit = distunit, file = trim(BndDistFileName), action = 'write', &
        &  form = 'unformatted', access = 'direct', &
        &  recl = rl)

      do iElem=1,temData%qVal(1)%nVals
        do iDir=1,qQQQ
          qVal(iDir) = temData%qVal(iDir)%val(iElem)
        enddo
        write(distunit, rec=iElem) qVal
      end do

      close(distUnit)
    endif ! temData%qVal(1)%nVals > 0
    ! ---------------------------------------------------------------------------

    ! ---------------------------------------------------------------------------
    ! Dump store_normal info to normals.lua
    anyNormals: if (temData%bc_normal(1)%nVals > 0) then
      tem_global%nProperties = tem_global%nProperties + 1
      property(tem_global%nProperties)%bitpos = prp_hasNormal
      property(tem_global%nProperties)%label  = 'has normal'
      property(tem_global%nProperties)%nElems = temData%bc_normal(1)%nVals
      call aot_out_open( conf, trim(confHeader%folder)//'normals.lua' )

      ! Dump logical variables that indicate whether the corresponding BC
      ! has normals stored for it
      call aot_out_open_table( conf, 'hasNormal' )
      do iBC = 1, actBC%nVals
        iUnibc = actBC%val(iBC)
        call aot_out_val( conf, &
          &               val = geometry%attribute%bc_uni_storenormal%val(iUnibc) )
      end do
      call aot_out_close_table(conf)
      call aot_out_close(conf)
      call destroy(actBC)
      ! -------------------------------------------------------------------------

      ! -------------------------------------------------------------------------
      ! Dump boundary normals to normals.lsb
      distUnit = newunit()
      BndDistFileName = trim(confHeader%folder) // 'normals' // EndianSuffix

      inquire(iolength = rl) bc_normal(:)
      open(unit = distunit, file = trim(BndDistFileName), action = 'write', &
        &  form = 'unformatted', access = 'direct', &
        &  recl = rl)

      do iElem=1,temData%bc_normal(1)%nVals
        do iDir=1,3
          bc_normal(iDir) = temData%bc_normal(iDir)%val(iElem)
        end do
        write(distunit, rec=iElem) bc_normal
      end do

      close(distUnit)
    end if anyNormals
    ! ---------------------------------------------------------------------------

    ! ---------------------------------------------------------------------------
    ! JQ: Dump property prp_solid into header file
    if ( temData%nSolids > 0 ) then
      tem_global%nProperties = tem_global%nProperties + 1
      ! JQ: add prp_solid property information
      ! propery 3 is has qval
      property(tem_global%nProperties)%bitpos = prp_solid
      property(tem_global%nProperties)%label  = 'solid'
      property(tem_global%nProperties)%nElems = temData%nSolids
    end if
    ! ---------------------------------------------------------------------------


    ! ---------------------------------------------------------------------------
    ! Dump color info to colors.lua
    if (temData%colors%nVals > 0) then
      tem_global%nProperties = tem_global%nProperties + 1

      property(tem_global%nProperties)%bitpos = prp_isColored
      property(tem_global%nProperties)%label  = 'is colored'
      property(tem_global%nProperties)%nElems = temData%colors%nVals

      call aot_out_open( conf, trim(confHeader%folder)//'colors.lua' )

      call aot_out_open_table( conf, 'color_label' )
      nColors = 0
      do iColor = 1,geometry%attribute%uni_name(sdr_Seed_object)%nVals
        col_name = geometry%attribute%uni_name(sdr_Seed_object)%val(iColor)
        if (trim(col_name) /= 'none') then
          write(logunit(2),*) 'Color ', trim(col_name), ' with a volume of ', &
            &                 temData%color_volume(iColor)
          write(logunit(2),*) 'in the mesh.'
          call aot_out_val( conf, val = trim(col_name) )
          nColors = nColors + 1
        end if
      end do
      call aot_out_close_table(conf)

      call aot_out_val( put_conf = conf, vname = 'nColors', &
        &               val = nColors )

      call aot_out_open_table( conf, 'color_fill' )
      do iColor = 1,geometry%attribute%uni_name(sdr_Seed_object)%nVals
        col_name = geometry%attribute%uni_name(sdr_Seed_object)%val(iColor)
        if (trim(col_name) /= 'none') then
          vpos = confHeader%subresolution%value_pos(iColor)
          cval = confHeader%subresolution%color_values%fill%val(vpos)
          call aot_out_val( conf, val = cval )
        end if
      end do
      call aot_out_close_table(conf)

      call aot_out_open_table( conf, 'color_void' )
      do iColor = 1,geometry%attribute%uni_name(sdr_Seed_object)%nVals
        col_name = geometry%attribute%uni_name(sdr_Seed_object)%val(iColor)
        if (trim(col_name) /= 'none') then
          vpos = confHeader%subresolution%value_pos(iColor)
          cval = confHeader%subresolution%color_values%void%val(vpos)
          call aot_out_val( conf, val = cval )
        end if
      end do
      call aot_out_close_table(conf)

      call aot_out_open_table( conf, 'color_volume' )
      do iColor = 1,geometry%attribute%uni_name(sdr_Seed_object)%nVals
        col_name = geometry%attribute%uni_name(sdr_Seed_object)%val(iColor)
        if (trim(col_name) /= 'none') then
          vpos = confHeader%subresolution%value_pos(iColor)
          cval = confHeader%subresolution%color_values%void%val(vpos)
          call aot_out_val( conf, val = temData%color_volume(iColor) )
        end if
      end do
      call aot_out_close_table(conf)


      call aot_out_close(conf)

      ! -------------------------------------------------------------------------
      ! Dump element colors to colors.ascii
      colunit = newunit()
      ColFileName = trim(confHeader%folder) // 'colors.ascii'

      inquire(iolength = rl) temData%colors%val(:,1)
      open( unit = colunit, file = trim(ColFileName), action = 'write', &
        &   form = 'unformatted', access = 'stream'                     )

      write(colunit) temData%colors%val

      close(colunit)
    end if
    ! ---------------------------------------------------------------------------


    ! ---------------------------------------------------------------------------
    ! Dump subresolution info to subresolution.lua
    if (temData%subres%nVals > 0) then
      tem_global%nProperties = tem_global%nProperties + 1

      property(tem_global%nProperties)%bitpos = prp_hasPolynomial
      property(tem_global%nProperties)%label  = 'has polynomial'
      property(tem_global%nProperties)%nElems = temData%subres%nVals

      call aot_out_open( conf, trim(confHeader%folder)//'subresolution.lua' )

      call aot_out_val( put_conf = conf,                               &
        &               vname    = 'polydegree',                       &
        &               val      = confHeader%subresolution%polydegree )
      call aot_out_val( put_conf = conf,                              &
        &               vname    = 'polyspace',                       &
        &               val      = confHeader%subresolution%polyspace )

      call aot_out_open_table( conf, 'projection' )
      call ply_prj_header_out( me = confheader%subresolution%project_header, &
        &                      conf = conf                                   )
      call aot_out_close_table(conf)

      call aot_out_open_table( conf, 'color_subres_nelems' )
      do iColor = 1,geometry%attribute%uni_name(sdr_Seed_object)%nVals
        col_name = geometry%attribute%uni_name(sdr_Seed_object)%val(iColor)
        if (trim(col_name) /= 'none') then
          call aot_out_val( conf, val = temData%color_subres_count(iColor) )
        end if
        ! Close the data files now, they were already written in create_target.
        if (temData%color_subres_count(iColor) > 0) then
          close( temData%color_subres_unit(iColor) )
        end if
      end do
      call aot_out_close_table(conf)

      call aot_out_close(conf)

      ! -------------------------------------------------------------------------
      ! Dump subresolved element colors to subres.ascii
      colunit = newunit()
      ColFileName = trim(confHeader%folder) // 'subres.ascii'

      inquire(iolength = rl) temData%subres%val(:,1)
      open( unit = colunit, file = trim(ColFileName), action = 'write', &
        &   form = 'unformatted', access = 'stream'                     )

      write(colunit) temData%subres%val

      close(colunit)
    end if
    ! ---------------------------------------------------------------------------

    if (temData%bc_ID(1)%nVals > 0) then
      call destroy(actBC)
    end if

    !mesh header info
    tem_global%BoundingCubeLength = geometry%universe%extent
    tem_global%Origin = geometry%universe%origin
    tem_global%effOrigin = temData%meshUniverse%global%effOrigin
    tem_global%effLength = temData%meshUniverse%global%effLength &
      &                  - temData%meshUniverse%global%effOrigin
    tem_global%nElems = temData%treeID%nVals
    tem_global%nParts = 1
    tem_global%myPart = 0
    tem_global%minLevel = temData%minlevel
    tem_global%maxLevel = temData%maxlevel
    tem_global%label = 'fluidTree'
    tem_global%comment = trim(confHeader%comment)
    tem_global%dirname = trim(confHeader%folder)

    ! allocate property list and propHeader infos
    allocate( tem_global%property( tem_global%nProperties ) )
    do iProp = 1, tem_global%nProperties
      tem_global%property(iProp) = property(iProp)
    end do
    deallocate( property )

    ! Dump global mesh info to header.lua.
    call dump_tem_global(tem_global)

    ! ---------------------------------------------------------------------------
    ! Dump mesh to elemlist.lsb or elemlist.msb
    meshunit = newunit()
    ElemFileName = trim(tem_global%dirname) // 'elemlist' // EndianSuffix

    inquire(iolength = rl) 0_long_k, 0_long_k
    open(unit = meshunit, file = trim(ElemFileName), action = 'write', &
      &  form = 'unformatted', access = 'direct', &
      &  recl = rl)

    do iElem = 1, temData%treeID%nVals
      write(meshunit, rec=iElem) temData%treeID%val(iElem), &
        & temData%PropertyBits%val(iElem)
    end do

    close(meshunit)

    write(logunit(1),*)'Seeder created mesh successfully!'

    call tem_stopTimer( timerHandle = timer_handle_dumping_mesh )

  end subroutine sdr_dump_treelm


end module sdr_proto2treelm_module
