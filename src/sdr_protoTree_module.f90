! Copyright (c) 2012-2016, 2020 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2012-2015, 2019 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2012 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2012, 2014, 2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2015 Samuel Ziegenberg
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
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
! ****************************************************************************
!> This module provides the description of the evolving tree, as it is built
!! up to match the given geometry.
!!
!! The concept used here, is a proposal for an algorithm, following a
!! level-wise iterative approach instead of an recursive one.
!! This is better suited for parallel processing.
!! Also it tries to follow more closely the principals from tem_construction
!! and the solver layout.
!! That is, treating seeder more like a solver in itself.
!!
!! For flooding a concept of wet sides on coarser levels is used, all elements
!! need to track just a single neighbor in each direction in this case.
!! And the coarse ghost elements take care of flooding their children
!! accordingly themselves.
module sdr_protoTree_module
  use, intrinsic :: iso_c_binding

  use env_module,            only: globalMaxLevels, long_k, rk, newunit,       &
    &                              isLittleEndian, pathLen, labelLen
  use tem_param_module,      only: qOffset, qQQQ, qInvDir, childCoord
  use treelmesh_module,      only: treelmesh_type
  use tem_debug_module,      only: main_debug
  use tem_global_module,     only: tem_global_type, dump_tem_global
  use tem_dyn_array_module,  only: PositionOfVal
  use tem_grow_array_module, only: grw_longArray_type, init, append, destroy,  &
    &                              truncate, grw_intArray_type
  use tem_restart_module,    only: tem_restart_type, tem_restart_writeHeader
  use tem_tools_module,      only: tem_horizontalSpacer
  use tem_logging_module,    only: logunit
  use tem_topology_module,   only: tem_directChildren, tem_FirstIdAtLevel,     &
    &                              tem_parentOf, tem_CoordOfId, tem_IdOfCoord
  use tem_varSys_module,     only: tem_varSys_type, tem_varSys_init, &
    &                              tem_varSys_append_stateVar, &
    &                              tem_varSys_proc_point, &
    &                              tem_varSys_proc_element
  use tem_varMap_module,     only: tem_create_varMap   
  use tem_time_module,       only: tem_time_type
  use tem_aux_module,        only: tem_abort
  use tem_solveHead_module,  only: tem_solverTag

  use sdr_geometry_module,   only: sdr_geometry_type, is_intersecting,         &
    &                              sdr_append_distanceRefineObject
  use sdr_spatialObj_module, only: periodicPlane
  use tem_cube_module,       only: tem_cube_type
  use sdr_attribute_module,  only: sdr_Boundary_Object, sdr_Seed_Object,       &
    &                              sdr_Refinement_Object,                      &
    &                              sdr_object_kinds_max,                       &
    &                              sdr_Fluidifyable_Object,                    &
    &                              sdr_noSolidification_Object,                &
    &                              sdr_any_bc_subresolution,                   &
    &                              sdr_any_bc_distanceRefine
  use sdr_node_module,       only: sdr_node_type, sdr_wetNeighborsFace,        &
    &                              intersectsBoundary_bit,  hasBoundary_bit,   &
    &                              isLeaf_bit, isTarget_bit, isFlooded_bit,    &
    &                              isFluidifyable_bit, isNoSolidification_bit, &
    &                              sdr_mark_floodNode, append, init, truncate, &
    &                              sdr_set_nodeProp_bit,                       &
    &                              sdr_clear_nodeProp_bit, sdr_nodeProp_btest, &
    &                              sdr_inHeritBnd_eligibleChildren,            &
    &                              sdr_append_childIntersectedObject,          &
    &                              sdr_intersectObjPos_type
  use sdr_config_module,     only: sdr_confHead_type
  use sdr_protoData_module,  only: sdr_append_protoVar

  use sdr_timer_module,      only: timer_handle_proto 
  use tem_timer_module,      only: tem_startTimer,tem_stopTimer 

  implicit none

  private

  public :: sdr_protoTree_type
  public :: sdr_build_protoTree
  public :: levelValues_type
  public :: sdr_write_proto_as_restart
  public :: sdr_neighbor_in_proto
  public :: sdr_node_neighbors


  !> The protoTree is used to describe the preliminary tree, before it
  !! is actually extended to the full information of the mesh.
  !!
  !! This extension will only be done for the actual leaf nodes, in the
  !! computational domain.
  type sdr_protoTree_type

    !> Keep track of the deepest level in the tree.
    integer :: nLevels

    !> Number of leaf nodes.
    integer :: nLeafNodes

    !> Number of flooded leaves.
    integer :: nFloodedLeaves

    !> List of all nodes in the tree, each node contains a list of the
    !! objects it is intersecting and a property bit mask.
    !! It is identified by its treeID.
    !!
    !! All components of Node are growing arrays, thus the corresponding
    !! component of node inode has to be accessed with component%val(inode).
    type(sdr_node_type) :: node

    !> The index of the first node on a given level in node%treeID%sorted.
    integer :: levelNode_first(0:globalMaxLevels)

    !> The index of the last node on a given level in node%treeID%sorted.
    integer :: levelNode_last(0:globalMaxLevels)

    !> Temporary array intersected objects of 8 children which will later
    !! be copied to intersected_object in node_type.
    !!
    !! This array is initialized in build_protoTree and destroyed
    !! after refine_leaf
    type(grw_intArray_type) :: child_intersected_object

  end type sdr_protoTree_type



  !> Auxilary data type to provide data on the current level iteration.
  type levelValues_type
    !> First treeID on current level
    integer(kind=long_k) :: ID_offset

    !> Node cube length on current level
    real(kind=rk) :: dx

    !> Number of the level itself.
    integer :: level
  end type levelValues_type


contains


  ! ****************************************************************************
  !> This routine builds the preliminary tree with geometry intersection and
  !! neighbor identification
  subroutine sdr_build_protoTree(me, geometry, header)
    ! --------------------------------------------------------------------------!
    !> preliminary tree created by this routine
    type(sdr_protoTree_type), intent(out) :: me
    !> type which contains all geometry object infos
    type(sdr_geometry_type), intent(inout) :: geometry
    !> some global information on solver name and version
    type(sdr_confHead_type), intent(inout) :: header
    ! --------------------------------------------------------------------------!
    integer(kind=long_k) :: parent_ID_offset
    real(kind=rk) :: parent_dx
    integer :: parent_coord(4)
    integer :: iLevel
    integer :: iParent
    integer :: iSide
    integer :: firstParent
    integer :: lastParent
    integer :: nodePos
    integer :: nNodes
    integer :: nColors
    integer :: iColor
    integer :: nBoundaries
    integer :: nSublevels
    logical :: distanceRefine
    integer, allocatable :: nodeprops(:)
    integer :: linkpos(6)
    logical :: nullify_parent_sl
    logical :: neigh_noLeaf, neigh_noSub, neigh_intersected
    logical :: testAll

    type(levelValues_type) :: leVal
    ! --------------------------------------------------------------------------!

    call tem_startTimer( timerHandle = timer_handle_proto ) 

    call tem_horizontalSpacer(funit=logunit(1))
    write(logunit(1),*) 'Building preliminary tree ...'
    me%nLeafNodes = 0
    me%nFloodedLeaves = 0
    me%levelNode_first = 1
    me%levelNode_last  = 0
    me%nLevels = 0

    ! Initialize temporary growing array for childrens
    ! Destroy it after refine_leaf routine
    call init(me%child_intersected_object)

    nColors = geometry%attribute%uni_name(sdr_Seed_object)%nVals
    nBoundaries = geometry%attribute%uni_name(sdr_Boundary_object)%nVals

    ! Set the number of sublevels only, if there is at least one active bc
    ! that requires it. Otherwise we are going to ignore subelement resolution.
    if (sdr_any_bc_subresolution(geometry%attribute)) then
      nSublevels = header%subresolution%nLevels
    else
      nSublevels = 0
    end if

    ! Set if distance refine is defined for any boundary attribute.
    ! Consider only boundary with distance_refine%reach_level>0
    distanceRefine = sdr_any_bc_distanceRefine(geometry%attribute)

    ! Initialize the list of all nodes
    call init( me             = me%node,                         &
      &        nColors        = nColors,                         &
      &        nSublevels     = nSublevels,                      &
      &        distanceRefine = distanceRefine,                  &
      &        none_color_id  = geometry%attribute%none_color_id )

    allocate(nodeprops(me%node%propLength))

    ! Define the root node (complete universe)
    nodeprops = 0

    ! By definition, this node intersects all geometric objects.
    ! so test all spatial objects only in the root node
    testAll = .true.

    ! The root node. It has a treeID of 0, and never is a subelement.
    call append(me           = me%node,         &
      &         treeID       = 0_long_k,        &
      &         PropertyBits = nodeprops,       &
      &         sublevel     = -1,              &
      &         minLevel     = header%minLevel, &
      &         pos          = nodePos          )

    ! Checking if there are any boundary objects here. To include periodic
    ! walls, we need to check all attribute kinds and can not rely on the
    ! counted kindpos, as that excludes the periodic objects.
    if (any(geometry%attribute%dynArray &
      &                       %val(:geometry%attribute%dynArray &
      &                                               %nVals)%kind &
      &     == sdr_Boundary_object)) then
      ! There are boundaries defined, as by definition, all boundaries are in
      ! the root node, the root node intersects boundaries, set the according
      ! bit to create its children subsequently.
      call sdr_set_nodeProp_bit(node  = me%node,              &
        &                       iNode = nodePos,              &
        &                       bit   = IntersectsBoundary_bit)

      ! No neighbors for this virtual node.
      do iSide=1,6
        me%node%linkpos(iSide)%val(nodePos) = 0
      end do

    else

      ! No boundaries defined at all, there is only the root element.
      ! (Which is also to be flooded unconditionally with all colors.)
      call sdr_set_nodeProp_bit(node  = me%node,  &
        &                       iNode = nodePos,  &
        &                       bit   = isLeaf_bit)
      do iColor=1,nColors
        call sdr_mark_floodNode(node           = me%node,           &
          &                     iNode          = nodePos,           &
          &                     nFloodedLeaves = me%nFloodedLeaves, &
          &                     color          = iColor             )
      end do

      ! The universe is periodic, and therefore all neighbors refer back to the
      ! root node itself.
      do iSide=1,6
        me%node%linkpos(iSide)%val(nodePos) = nodePos
      end do
      ! Increase the counter of leaf nodes
      me%nLeafNodes = me%nLeafNodes + 1

    end if

    ! Mark the boundaries of level 0.
    me%levelNode_first(0) = nodePos
    me%levelNode_last(0) = nodePos

    ! Set the length of the level cube length to the complete bounding cube
    ! length initially.
    leVal%dx = geometry%universe%extent
    leVal%ID_offset = 0_long_k

    levelLoop: do iLevel = 1, globalMaxlevels
      write(logunit(2),*) 'current Level ', iLevel
      ! The first node on the current level has to be one after all the nodes
      ! already in the mesh.
      me%levelNode_first(iLevel:) = me%node%nNodes+1

      firstParent = me%levelNode_first(iLevel-1)
      lastParent = me%levelNode_last(iLevel-1)
      parent_ID_offset = leVal%ID_offset ! First treeID on parent level
      parent_dx = leVal%dx
      ! (save from previous iteration)

      ! Set some auxilary data describing the current level.
      leVal%dx = 0.5_rk * leVal%dx ! length of cube nodes on this level
      leVal%ID_offset = tem_FirstIdAtLevel(iLevel) ! first treeID on this level
      leVal%level = iLevel ! level count

      ! Set hasBoundary properties for neighbors of intersected 
      ! leaf bit on parent level to inHerit to children before looping
      ! over parent to create children. 
      ! Also create sphere object if distance function
      ! is defined for this boundary and append sphere object
      ! to list of spatial objects.
      do iParent = firstParent, lastParent
        ! Coordinates of parent
        parent_coord = tem_CoordOfId(                                &
          &                    treeID = me%node%treeID%val(iParent), &
          &                    offset = parent_ID_offset             )

        if ( sdr_nodeProp_btest(node  = me%node, &
          &                     iNode = iParent, &
          &                     bit   = intersectsBoundary_bit) ) then

          ! if intersected boundary and leaf node then set
          ! neighbors of this node to hasBoundary bit
          if ( sdr_nodeProp_btest(node  = me%node,  &
            &                     iNode = iParent,  &
            &                     bit   = isLeaf_bit) ) then
            call sdr_mark_neighborHasBnd(proto = me,         &
              &                          coord = parent_coord)
          end if

          ! Create sphere object and attribute for distance function and 
          ! append to list of spatial objects if parent in intersected boundary
          call sdr_append_distanceRefineObject(                            &
            &             coord              = parent_coord,               &
            &             dx                 = parent_dx,                  &
            &             iLevel             = iLevel-1,                   &
            &             geometry           = geometry,                   & 
            &             intersected_first  = me%node%userObjPos          &
            &                                         %val(iParent)%first, &
            &             intersected_last   = me%node%userObjPos          &
            &                                         %val(iParent)%last,  &
            &             intersected_object = me%node%intersected_object  )
        end if ! parent is intersectedBoundary   
      end do ! parent loop

      ! Iterate over all elements on the previous level.
      parElemLoop: do iParent = firstParent, lastParent

        if ( sdr_nodeProp_btest(node  = me%node,  &
          &                     iNode = iParent,  &
          &                     bit   = isLeaf_bit) ) then

          ! At this point the domain boundaries have been resolved down to the
          ! required level in this location of the mesh.

          ! For flooding we need the neighbors of all leaf-nodes, and children
          ! of non-leaves.
          ! We can exploit here, that the 8 children are always placed
          ! consecutively one after the other in the node list an thus, storing
          ! the position of the first child is sufficient to find them all.
          ! 6 integers to store the positions of linked elements are therefore
          ! sufficient, as neighbors are not needed for virtual nodes.
          ! Reached a leaf node in the prototree, lookup it's neighbors.
          ! This is possible here, as prototree neighbors are either on the
          ! current level, or on a previous level.
          linkpos = sdr_node_neighbors( me       = me,               &
            &                       level_offset = parent_ID_offset, &
            &                       iNode        = iParent           )
          do iSide=1,6
            me%node%linkpos(iSide)%val(iParent) = linkpos(iSide)
          end do

          ! If the current node is flooded, wet all the neighboring sides
          ! by setting the appropiate bits in the neighbors properties.
          ! This avoids special treatment in the flooding routine for the
          ! initially flooded elements.
          call sdr_wetNeighborsFace(node  = me%node, &
            &                       iNode = iParent  )

        else

          ! The parent is not a leaf node, will add all eight children at the
          ! end of the current node list, thus the first one will be at position
          ! nNodes + 1.
          ! Store this vertical link here, as it will later on be required
          ! during flooding.
          me%node%linkpos(1)%val(iParent) = me%node%nNodes + 1

          ! Flag to indicate, whether the parent sublevel needs to be restored
          ! to a non-negative value, after children where created.
          nullify_parent_sl = .false.

          if ( sdr_nodeProp_btest(node  = me%node,    &
            &                     iNode = iParent,    &
            &                     bit   = isTarget_bit) &
            &  .and. geometry%smoothbounds              ) then
            ! If the parent was a target, we need to check the neighbors for
            ! their status, to ensure proper resolution, in the proximity of
            ! highly resolved non-subelement boundaries.
            ! Obviously this can only happen, if subresolution is active, and
            ! therefore the sublevel array is available.
            ! By setting smoothbounds to false this can be deactivated, however
            ! in this case it might happen, that the subresolved element
            ! neighbors to multiple boundary conditions on its sides and it is
            ! not quite clear which ones to use.
            linkpos = sdr_node_neighbors( me           = me,           &
              &                       level_offset = parent_ID_offset, &
              &                       iNode        = iParent           )
            do iSide = 1,6
              ! If any adjacent element is not a leaf but intersected by
              ! boundaries, it will be further refined. In this case we also
              ! need to refine the current parent further, that is, we need to
              ! remove its target bit, and turn it into a non-subelement node
              ! for the children (sublevel < 0).
              ! However, this should only be done, if the intersected neighbor
              ! itself is not a subelement. To avoid neighboring target elements
              ! to cause infinite refinement, we therefore need to detect
              ! neighboring target elements as well. As the sublevel has no
              ! relevance after creating the chilren, we will temporaly set the
              ! sublevel of the current parent to a negative value and restore
              ! a positive one afterwards again (for future neighbor checks).
              ! Though, this might be a little confusing, the only other option
              ! would be to introduce yet another bit to set for this case.
              neigh_noLeaf = .not. sdr_nodeProp_btest( node  = me%node,        &
                &                                      iNode = linkpos(iSide), &
                &                                      bit   = isLeaf_bit      )
              neigh_intersected = sdr_nodeProp_btest(              &
                &                   node  = me%node,               &
                &                   iNode = linkpos(iSide),        &
                &                   bit   = intersectsBoundary_bit )
              neigh_noSub = me%node%sublevel%val(linkpos(iSide)) < 0

              if (neigh_noSub .and. neigh_intersected .and. neigh_noLeaf) then
                ! Found a neighbor to the current parent, that will be further
                ! refined due to boundaries. Remove the target bit from it and
                ! pretend a negative sublevel for the chilren.
                nullify_parent_sl = .true.
                me%node%sublevel%val(iParent) = -1
                call sdr_clear_nodeProp_bit( node  = me%node,     &
                  &                          iNode = iParent,     &
                  &                          bit   = isTarget_bit )
                ! We can end the loop as soon, as any neighbor with these
                ! properties was found...
                EXIT
              end if

            end do
          end if

          ! Add all eight children.
          call create_children(me       = me,             &
            &                  parent   = iParent,        &
            &                  geometry = geometry,       &
            &                  leVal    = leVal,          &
            &                  testAll  = testAll,        &
            &                  minlevel = header%minlevel )
          
          if (nullify_parent_sl) then
            ! If we removed the target bit from our parent above, we need to set
            ! its sublevel to a non-negative value here again.
            ! Otherwise neighboring target nodes would not recognize this node
            ! as one, which is only refined due to neighbors.
            me%node%sublevel%val(iParent) = 0
            ! Reset the parent to a common virtual node again, thus it should
            ! have a 0 minbcid again.
            me%node%minBCID%val(iParent) = 0
          end if

        end if ! parent is leaf

      end do parElemLoop


      me%levelNode_last(iLevel:) = me%node%nNodes

      ! The number of elements on the current level, compute this count with
      ! the difference of the first and last index of nodes on the current
      ! level.
      nNodes = me%levelNode_last(iLevel) - me%levelNode_first(iLevel) + 1
      write(logunit(3),*) '  * nodes on this level: ', nNodes
      write(logunit(3),*) ''

      ! Leave the loop, if no new nodes were created in this iteration.
      if (nNodes == 0) exit

      ! Update the number of levels in the tree, as the current level contains
      ! new elements.
      me%nLevels = iLevel

      ! after root node deactivate testAll
      testAll = .false.

    end do levelLoop

    ! Found all leaves and resolved boundaries to the requested resolution.
    ! Free memory of growing arrays, if they are longer than the actual
    ! data.
    call truncate(me%node)

    write(logunit(1),*) '                          ... done with level: ', &
      &                 me%nLevels
    write(logunit(1),*) 'Number of leaves:      ', me%nLeafNodes
    write(logunit(1),*) 'Total number of nodes: ', me%node%nNodes

    write(logunit(2),*) 'Number of flooded leaves while building tree: ', &
      &                 me%nFloodedLeaves

    write(logunit(1),*) 'Memory unused in parent intersected_object: ', &
      &                 me%node%memLeft_userObj
    if ( me%nFloodedLeaves == 0 ) then
      write(logunit(0),*) 'Error: There are no flooded leaves. '
      write(logunit(0),*) &
        &  'PLACE SEED AT NON-SOLID ELEMENT. Terminating seeder...'
      call tem_abort()
    end if
    write(logunit(1),*) ''

    ! Finished the prototree with all necessary data for flooding.
    ! We can now go on and identify the elements belonging to the domain with
    ! the flooding algorithm.

    call tem_stopTimer( timerHandle = timer_handle_proto )

  end subroutine sdr_build_protoTree
  ! ****************************************************************************


  ! ****************************************************************************
  !> This routine creates children for each parent if children
  !! intersect with boundary object.
  !!
  !! First loop over 8 children, and test for intersection of each child cube
  !! with the geometry objects inherited from the parent node.
  !! Then check for various object kinds, that might be intersected:
  !! If BOUNDARY objects are intersected, record the minimal bcID for later,
  !! boundaries are only marked as a leaf if the maxlevel has been reached.
  !! If NO BOUNDARY are intersected, the refinement can stop early here, and
  !! the node is marked as leaf. This avoids overly many elements before
  !! flooding.
  !! If a leaf is intersecting a SEED object, mark it already as flooded here.
  !! Children that do not intersect any objects do not need to be refined
  !! further at this stage and are marked as leaf nodes.
  subroutine create_children(me, parent, geometry, leval, testAll, minlevel)
    ! --------------------------------------------------------------------------!
    !> preliminary tree on which childern are created
    type(sdr_protoTree_type), intent(inout) :: me
    !> Position of parent node on the growing array of node_treeID and node_data
    !! in preliminary tree
    integer, intent(in) :: parent
    !> type which contains all geometry object infos
    type(sdr_geometry_type), intent(in) :: geometry
    !> contains information on current level on which children are created
    type(levelValues_type), intent(in) :: leVal
    !> testAll objects only for root node
    logical, intent(in) :: testAll
    integer, intent(in) :: minlevel !minimum refinement for fluids
    ! --------------------------------------------------------------------------!
    integer :: iChild, iObject
    integer :: nObjects ! number of objects need to intersect with
    integer(kind=long_k) :: treeID(8) ! treeID of eight children
    integer :: node_coord(4)
    integer :: firstchild_coord(4)
    integer :: node_pos ! position in node list ( me%node )
    integer :: obj_pos
    integer :: attr_pos
    integer :: attr_kind
    integer :: maxLevel ! the ultimate level this child should resolve
    integer :: min_BCID
    integer :: n_ni_seeds
    integer :: seed_color(me%node%nColors)
    integer :: propbits(me%node%propLength)
    integer :: iColor
    integer :: col_int, col_bit
    ! number of objects in each kind for one child
    integer :: objkind_count(sdr_object_kinds_max)
    integer :: bc_color
    integer :: bcid
    integer :: child_sublevel
    logical :: reached_leaf
    logical :: reached_maxlevel
    logical :: subresolution
    logical :: color_seeded(me%node%nColors)
    logical :: color_intersected(me%node%nColors)
    logical :: nonecolor_intersected
    logical :: child_hasBnd(8)
    type(sdr_intersectObjPos_type) :: child_objPos(8)
    integer :: intersected_first, intersected_last
    integer :: child_nodePos(8)
    integer :: child_nObjects(8)
    integer :: memLeft

    type(tem_cube_type) :: node_cube
    ! --------------------------------------------------------------------------!

    ! If parent has hasBoundary_bit then inherit this property to its children
    child_hasBnd = sdr_inHeritBnd_eligibleChildren(me%node, parent)

    ! Get treeID of eight children
    treeID(1:8) = tem_directChildren(me%node%treeID%val(parent))

    if (me%node%subelement_resolution > 0) then
      child_sublevel = me%node%sublevel%val(parent) - 1
    else
      ! If there is no subelement resolution the sublevels are not available.
      ! However, we store a negative value here in child_sublevel to avoid
      ! further checks for subelement_resolution further down.
      child_sublevel = -1
    end if

    ! Get the length of children cube. all children cubes have same length.
    node_cube%extent = leVal%dx
    node_cube%halfwidth = 0.5_rk * leVal%dx

    if (testAll) then
      nObjects = geometry%nUserObjs
      intersected_first = 1
      intersected_last = nObjects
    else
      intersected_first = me%node%userObjPos%val(parent)%first
      intersected_last = me%node%userObjPos%val(parent)%last
      nObjects = intersected_last - intersected_first + 1
    end if

    ! Coordinate of childs can be computed with an offset from
    ! first child using childCoord parameter
    firstchild_coord = tem_coordOfId(treeID = treeID(1),     &
      &                              offset = leVal%ID_offset)
 
    ! Initialize counter for child_intersected_first and last     
    child_objPos(:)%first = 1
    child_objPos(:)%last = 0
    me%child_intersected_object%nVals = 0

    childLoop: do iChild = 1, 8
      !>@todo HK: The content of the childLoop should probably move into its
      !!          own subroutine, it could then also be used to define the root
      !!          node (treeID=0). Though this might just add a 9th check on
      !!          all defined objects in most cases. Needs further thinking.
      ! Set the flag on the maximal level in the intersected boundary objects
      ! initially to 0.
      maxLevel = minlevel
      ! Set the counts for all object kinds to 0.
      objkind_count = 0

      ! the first intersected object of this child has to be one after
      ! all the previous childs intersected objects
      child_objPos(iChild:)%first = me%child_intersected_object%nVals + 1

      ! Define the children cube to intersect with.
      node_coord(1:3) = firstchild_coord(1:3) + childCoord(iChild,:) 
      node_cube%origin = leVal%dx * node_coord(1:3) + geometry%universe%origin
      node_cube%endPnt = leVal%dx * (node_coord(1:3)+1) + geometry%universe%origin
      node_cube%center = node_cube%origin + node_cube%halfwidth

      ! Append this child to the tree.
      propBits = 0
      call append(me           = me%node,        &
        &         treeID       = treeID(iChild), &
        &         PropertyBits = propBits,       &
        &         sublevel     = child_sublevel, &
        &         minlevel     = minlevel,       &
        &         pos          = node_pos        )

      ! position of this node in node list
      child_nodePos(iChild) = node_pos  

      ! if this child hasBnd, set the hasBoundary_bit 
      ! and copy the directions which has boundary from parent.
      ! if parent has boundary in certain direction, its child should have 
      ! the save since hasBoundary is set by the leaf of intersected 
      ! boundary node
      if (child_hasBnd(iChild)) then
        call sdr_set_nodeProp_bit(node  = me%node,        &
          &                       iNode = node_pos,       &
          &                       bit   = hasBoundary_bit )
        me%node%hasBndBits%val(node_pos) = me%node%hasBndBits%val(parent)  
      end if

      ! Trumping rule to set boundaryID for each direction when a node
      ! is intersected by more than one boundary object.
      min_bcID = huge(min_bcID)
      bcid = 0

      color_seeded = .false.
      color_intersected = .false.

      ! Check for subresolution of this node.
      ! Subresolution is only activated for non-periodic boundaries and those
      ! with a specific color. Boundaries that should apply to all colors can
      ! not be subresolved, as it is unclear what values should be used outside
      ! the computational domain.
      ! Subresolution becomes necessary if any intersected boundary requests it.
      ! If a boundary should be subresolved or not is given by its attribute.
      subresolution = .false.
      nonecolor_intersected = .false.

      ! Loop over the objects which intersected with parent
      objLoop: do iObject = intersected_first, intersected_last
        if (testAll) then
          obj_pos = iObject
        else
          obj_pos = me%node%intersected_object%val(iObject)
        end if

        ! Test for intersection of geometric object and node_cube.
        if (is_intersecting(node_cube, geometry, obj_pos)) then
          ! If there is an intersection, add the object to the
          ! intersected_object list of the child.
          call append( me%child_intersected_object, obj_pos )
          attr_pos = geometry%spatialObj%val(obj_pos)%attribute_position
          attr_kind = geometry%attribute%dynArray%val(attr_pos)%kind
          objkind_count(attr_kind) = objkind_count(attr_kind) + 1

          select case(geometry%attribute%dynArray%val(attr_pos)%kind)
          case (sdr_Boundary_object)
            bcid = geometry%attribute%dynArray%val(attr_pos)%id
            min_bcID = min(min_bcID, bcid)
            ! Exclude periodic boundaries from the color check.
            if (bcid > 0) then
              bc_color = geometry%attribute%bc_color_id(bcid)
              if (bc_color < 0) then
                ! Negative values are used for boundaries, which should
                ! affect all colors.
                color_intersected = .true.
                ! The nonecolor is treated somewhat special as, a node that is
                ! intersected by this kind of boundary should not be
                ! subresolved.
                nonecolor_intersected = .true.
              else if (bc_color > 0) then
                ! This boundary should be applied to a specific color.
                color_intersected(bc_color) = .true.
                ! If there is any boundary with a specific color and
                ! required subresolution, this element needs to be resolved
                ! beyond the levels of the mesh.
                subresolution = (subresolution                             &
                  &              .or. geometry%attribute%dynArray          &
                  &                           %val(attr_pos)%subresolution )
              end if
            else
              ! Periodic boundaries always affect all colors.
              color_intersected = .true.
            end if

          case (sdr_Seed_object)
            color_seeded(geometry%attribute%dynArray%val(attr_pos)%uni_id) &
              &  = .true.
          end select

          ! Update the maximal level of intersected objects.
          ! ( Always refine objects down to maximal level of all objects
          ! intersecting the cube. That is, refinementboxes might enforce
          ! a higher boundary resolution than the boundary object itself.
          ! Objects, that should have no influence on the resolution, like
          ! seeds, should have set their level low enough e.g. 0).
          ! KM: Do not refine periodic boundaries, Refining periodic boundaries
          ! causes issue in finding periodic neighbor between coarser fluid
          ! node and refined periodic boundary node.
          if (bcid /= -1) then
            maxLevel = max(maxlevel, &
              &            geometry%attribute%dynArray%val(attr_pos)%level)
          end if  
        end if ! is_intersecting
      end do objLoop

      ! ------------------------------------------------------------------------!
      ! Set child intersected object first and last in child_intersected_object
      child_objPos(iChild:)%last = me%child_intersected_object%nVals
      child_nObjects(iChild) = child_objPos(iChild)%last    &
        &                    - child_objPos(iChild)%first + 1 
      ! ------------------------------------------------------------------------!

      ! Check whether this child intersects any boundary object.
      ! Update intersectsboundary_bit if necessary.
      ! For the proto tree only boundaries need to be resolved fully,
      ! stop everything else as early as possible.
      ! Thus, reached leaf depends on the maxlevel only if there are
      ! boundaries intersected it is true in all other cases.
      if (objkind_count(sdr_Boundary_object) > 0 ) then
        call sdr_set_nodeProp_bit( node  = me%node,               &
          &                        iNode = node_pos,              &
          &                        bit   = intersectsBoundary_bit )
        reached_maxlevel = (leval%level == globalmaxlevels)
        reached_leaf = ( (maxLevel <= leval%level) .or. reached_maxlevel)

        ! Only do a subresolution, if the element is not intersected by a BC
        ! with the 'none' color. BCs with the 'none' color restrict the overall
        ! computational domain, and elements intersected by it, will not be
        ! part of the domain. Therefore, there is no need for subresolution
        ! here.
        subresolution = ( subresolution .and. (.not. nonecolor_intersected) &
          &               .and. (me%node%subelement_resolution > 0)         )

        leaf: if (reached_leaf) then
          sub: if (child_sublevel < 0) then
            ! Hit either a leaf or target element (node is not subelement yet).
            me%node%minBCID%val(node_pos) = min_bcID

            if ( subresolution .and. (.not. reached_maxlevel) &
              &                .and. min_bcID >= 0            ) then
              ! There is a boundary that needs to be refined on a subelement
              ! resolution, this can only be true, if subelement_resolution is
              ! greater than 0, and thus, the sublevel data is available.
              ! Mark this node as target, and set the sublevels according to the
              ! subelement_resolution.
              ! Subresolution is not possible if the maximal level was already
              ! reached, in this case, the element is treated like a normal
              ! leaf node.
              ! It is also not possible periodic boundaries, if the element
              ! intersects a periodic boundary, subresolution will be ignored.
              ! (Periodic boundaries trump all other boundaries)
              call sdr_set_nodeProp_bit( node  = me%node,     &
                &                        iNode = node_pos,    &
                &                        bit   = isTarget_bit )
              me%node%sublevel%val(node_pos) = me%node%subelement_resolution
              reached_leaf = .false.
            end if

          else sub

            ! The parent is a target element or already a subelement.
            ! Decide end of refinement independent of maxLevel and use the
            ! sublevel counter instead. Never refine beyond the maximal allowed
            ! total number of elements.
            if ( (child_sublevel > 0) .and. (.not. reached_maxlevel) ) then
              reached_leaf = .false.
            else
              ! Reached the required subelement resolution, finally this can
              ! be marked as a leaf node, also store the minbcid for this
              ! element now.
              me%node%minBCID%val(node_pos) = min_bcID
            end if
          end if sub
        end if leaf

      else
        ! If there is no boundary intersected, no further refinement
        ! is necessary.
        reached_leaf = .true.

      end if

      ! Check whether this child intersects any noSolidification object.
      ! Update noSolidification_bit if necessary.
      if (objkind_count(sdr_noSolidification_object) > 0 ) then
        call sdr_set_nodeProp_bit(node  = me%node,              &
          &                       iNode = node_pos,             &
          &                       bit   = isNoSolidification_bit)
      end if

      ! Set count for non-intersected leaves to 0, to allow their counting
      ! in the loop.
      n_ni_seeds = 0

      ! Setting color specific properties according to the intersected objects
      ! now (only necessary if there is actually something intersected):
      if (child_nObjects(iChild) > 0) then

        ! Loop over colors to find color dependent properties of this node.
        do iColor=1,me%node%nColors
          if (color_intersected(iColor)) then

            ! The node is intersected by a boundary that is relevant for
            ! iColor. Set the according bit.
            col_int = (iColor-1) / me%node%bytes_per_int + 1
            col_bit = mod(iColor-1, me%node%bytes_per_int)*8
            me%node%PropertyBits%val(col_int,Node_pos) &
              &   = ibset(me%node%PropertyBits%val(col_int,Node_pos), col_bit)

          else if (color_seeded(iColor)) then

            ! Not intersected by relevant boundary, but seeded.
            ! Add to list of possible seeds. If this node is actually a leaf,
            ! it will be marked as flooded in that color.
            n_ni_seeds = n_ni_seeds + 1
            seed_color(n_ni_seeds) = iColor

          end if
        end do

      end if ! this node intersects with any object


      if (reached_Leaf) then
        ! This node is actually a leaf, mark it as such and flood it in
        ! seeded colors as needed.
        call mark_leafNode(me       = me,     &
          &                node_pos = node_pos)
        do iColor=1,n_ni_seeds
          ! For all non-intersected seeds, mark the leaf node as flooded.
          call sdr_mark_floodNode(node           = me%node,           &
            &                     iNode          = node_pos,          &
            &                     nFloodedLeaves = me%nFloodedLeaves, &
            &                     color          = seed_color(iColor) )
        end do
      end if

    end do childLoop


    ! Now copy temporary growing array of intersected objects of child
    ! to actual array of intersected_objects
    if (any(child_nObjects > 0)) then
      call sdr_append_childIntersectedObject(                          &
        &      node                     = me%node,                     &
        &      parent                   = parent,                      &
        &      testAll                  = testAll,                     &
        &      intersected_object       = me%node%intersected_object,  &
        &      grwObjPos                = me%node%userObjPos,          &
        &      child_nodePos            = child_nodePos,               &
        &      child_nObjects           = child_nObjects,              &
        &      child_intersected_object = me%child_intersected_object, &
        &      child_objPos             = child_objPos,                &
        &      memLeft                  = memLeft                      )

      ! track memory unused from parent intersected object list  
      me%node%memLeft_userObj = me%node%memLeft_userObj + memLeft  
    end if  

  end subroutine create_children
  ! ****************************************************************************


  ! ****************************************************************************
  !> Small helping routine to keep track of leaf nodes.
  !!
  !!@todo HK: Originally this operation was inlined in create_children in
  !!          a few places. Need to check for performance impact? Readability
  !!          was not hurt too much by this, but it was changed into a flag and
  !!          final decision, which in my opinion did not improve the
  !!          readability as the reader has to track down the location where the
  !!          final decision is taken, and introduces an additional if.
  !!          It might be that the compiler is smart enough on both variants
  !!          but anyway, I would prefer such a small routine if you think
  !!          inlining is too burdensome.
  !!          We should have an eye on the performance, though.
  subroutine mark_leafNode(me, node_pos)
    ! --------------------------------------------------------------------------!
    !> Property bit mask to set the leaf bit in
    type(sdr_protoTree_type), intent(inout) :: me  
    integer, intent(in) :: node_pos
    ! --------------------------------------------------------------------------!
    me%nLeafNodes = me%nLeafNodes + 1
    call sdr_set_nodeProp_bit(node  = me%node,  &
      &                       iNode = node_pos, &
      &                       bit   = isLeaf_bit)
  end subroutine mark_leafNode
  ! ****************************************************************************


  ! ****************************************************************************
  !> This routine identifies the 6 direct neighbors of a node in the prototree
  !!
  !!@todo We could cut the lookups of neighbors in the complete tree
  !!      by one half, if we set the siblings within the same
  !!      direct parent by their in-node relation (all 8 siblings within a
  !!      parent either exist or do not exist at all).
  !!      We than would have to iterate only over the remaining
  !!      three outer sides of each node here.
  function sdr_node_neighbors(me, level_offset, iNode, coord) result(neighbors)
    ! --------------------------------------------------------------------------!
    !> neighbors are identified for this tree and neighbor of each node
    !! are stored at link_pos of each node in node_data type
    type(sdr_protoTree_type), intent(in) :: me
    !> First treeID on level
    integer(kind=long_k), optional, intent(in) :: level_offset
    !> Node position in protoTree
    integer, intent(in) :: iNode
    !> if coord is present, no need to compute using tem_coordOfID
    integer, intent(in), optional :: coord(4)
    integer :: neighbors(6)
    ! --------------------------------------------------------------------------!
    integer :: iLevel
    integer :: iDir
    integer :: iSide
    integer :: iNeighbor
    integer :: offset(4)
    integer(kind=long_k) :: neighbor_tID
    integer :: neighbor_pos
    integer :: neighbor_level
    integer :: coord_loc(4)
    ! --------------------------------------------------------------------------!

    ! 3D offset to indicate neighbor to look up.
    offset = 0

    if (present(coord)) then
      coord_loc = coord
    else 
      coord_loc = tem_coordOfId(treeID = me%node%treeID%val(iNode), &
        &                       offset = level_offset)
    end if
    iLevel = coord_loc(4)

    dirLoop: do iDir=1,3

      sideLoop: do iSide=1,2

        offset(iDir) = iSide*2-3          ! offset = -1 or 1
        iNeighbor = (iDir-1)*2 + iSide    ! iNeighbor = 1, 2,..., 6
        neighbor_tID = tem_IdOfCoord(coord = coord_loc+offset, &
          &                          offset = level_offset)

        ! Find the neighbor on the current level or on any
        ! one above (move one level up, when the neighbor is not found).
        neighLoop: do neighbor_level=iLevel,0,-1
          neighbor_pos = PositionOfVal(                                & 
            &              me    = me%node%treeID,                     &
            &              val   = neighbor_tID,                       &
            &              lower = me%levelNode_first(neighbor_level), &
            &              upper = me%levelNode_last(neighbor_level)   )
          if (neighbor_pos > 0) then
            ! Found the neighbor, put it into the list of neighbors
            neighbors(iNeighbor) = neighbor_pos
            ! Leave the loop
            exit
          else
            neighbor_tID = tem_ParentOf(neighbor_tID)
          end if
        end do neighLoop

      end do sideLoop
      ! reset offset in current direction to 0.
      offset(iDir) = 0

    end do dirLoop

  end function sdr_node_neighbors
  ! ****************************************************************************


  ! ****************************************************************************
  !> This routine marks 26 direct neighbors as has boundary bit
  subroutine sdr_mark_neighborHasBnd(proto, coord)
    ! --------------------------------------------------------------------------!
    !> neighbors are identified for this tree and neighbors
    !! as marked with hasBoundary_bit
    type(sdr_protoTree_type), intent(inout) :: proto
    !> Coordinate of current node
    integer, intent(in) :: coord(4)
    ! --------------------------------------------------------------------------!
    integer :: iDir
    integer :: neighbor_level
    integer :: neighbor_pos
    ! --------------------------------------------------------------------------!

    ! loop over all 26 directions
    do iDir = 1, qQQQ
      ! get position of neighbor in the protoTree which might me
      ! from different level
      neighbor_pos = sdr_neighbor_in_proto( proto, coord, iDir,  &
        &                                   neighbor_level )
      
      call sdr_set_nodeProp_bit(node  = proto%node,    &
        &                       iNode = neighbor_pos,  &
        &                       bit   = hasBoundary_bit)


      ! From neighbor node the current intersected boundary node
      ! is at inverse direction of iDir
      proto%node%hasBndBits%val(neighbor_pos) =                &
        & ibset( proto%node%hasBndBits%val(neighbor_pos), qInvDir(iDir) )
    end do    
  end subroutine sdr_mark_neighborHasBnd
  ! ****************************************************************************


  ! ****************************************************************************
  !> Write current leaves of the prototree as treelm restart.
  !!
  !! This routine is mainly for debugging, to allow the visualization of the
  !! tree with harvester. The mesh will not contain any properties, instead
  !! additional data is provided as restart.
  !! Note that recursive output is used, to write each leaf node on its own,
  !! this output is thus really only useful for debugging, and might take very
  !! long for large meshes, however memory consumption should be fairly small.
  !! @todo SDR_write_proto_as_restart is a explicitly serial routine, deploying
  !!       its own writing to save memory. Therefore this is some code
  !!       duplication from treelmesh_module::dump_treelmesh.
  subroutine sdr_write_proto_as_restart(proto, geometry, itime, header, &
    &                                   prefix)
    ! --------------------------------------------------------------------------!
    !> The prototree to output.
    type(sdr_protoTree_type), intent(in) :: proto
    !> Bounding cube, the prototree lives in.
    type(sdr_geometry_type), intent(in) :: geometry
    !> wave number in character format to prepend to filenames and timestamp
    integer, intent(in) :: itime
    !> some global information on solver name and version
    type(sdr_confHead_type), intent(in) :: header
    !> prefix for filenames
    character(len=*), optional, intent(in) :: prefix
    ! --------------------------------------------------------------------------!
    type(treelmesh_type) :: dummy_treelmesh
    type(tem_global_type) :: tem_global
    type(tem_varsys_type) :: varsys
    type(tem_restart_type) :: restart
    type(tem_time_type) :: timing
    character(len=PathLen) :: ElemFileName
    character(len=4) :: EndianSuffix
    character(len=labelLen) :: charWave
    character(len=labelLen), allocatable :: varLabel(:)
    !character(len=labelLen) :: varLabel
    real(kind=rk) :: protodata(11+proto%node%nColors)
    integer :: meshunit
    integer :: restunit
    integer :: rl
    integer :: iLeaf
    integer :: iColor
    character(len=labelLen) :: prefix_loc
    procedure(tem_varSys_proc_point), pointer :: get_point => null()
    procedure(tem_varSys_proc_element), pointer :: get_element => null()
    ! --------------------------------------------------------------------------!

    if (isLittleEndian) then
      EndianSuffix = '.lsb'
    else
      EndianSuffix = '.msb'
    end if

    if(present(prefix)) then
      prefix_loc = trim(prefix)
    else
      prefix_loc = ''
    endif
    
    write(charwave,*) itime
    charwave=trim(adjustl(charwave))
    timing%iter = itime
    timing%sim = real(itime, kind=rk)

    tem_global%BoundingCubeLength = geometry%universe%extent
    tem_global%Origin = geometry%universe%origin
    tem_global%EffLength = tem_global%BoundingCubeLength
    tem_global%EffOrigin = tem_global%Origin
    tem_global%nElems = proto%nLeafNodes
    tem_global%nParts = 1
    tem_global%myPart = 0
    tem_global%minLevel = 0
    tem_global%maxLevel = proto%nLevels
    tem_global%label = 'ProtoTree'
    tem_global%comment = 'A ProtoTree from Seeder.'
    tem_global%nProperties = 0 ! No properties, everything in restart data.
    tem_global%dirname = trim(main_debug%debugMesh)//trim(prefix_loc)&
      &                 //trim(charWave) // '_mesh_'
    dummy_treelmesh%global = tem_global

    allocate( tem_global%property( tem_global%nProperties ) )

    ! Initialize varsys with a length of 16 (there will be at least 12 vars)
    call tem_varsys_init( me         = varsys,     &
      &                   systemName = 'meshInfo', &
      &                   length     = 16          )

    ! proto variables dumped as restart output
    allocate(varLabel(11+proto%node%nColors))
    varLabel(1:11) = (/ 'refLevel          ', 'isFlooded         ', &
      &                 'nObjects          ', 'nSeeds            ', &
      &                 'nRefines          ', 'nFluidifyable     ', &
      &                 'nBounds           ', 'trumping_bound    ', &
      &                 'target_level      ', 'bcAttr_pos_minBCID', &
      &                 'hasBnd            '                       /) 

    do iColor=1,proto%node%nColors
      ! 2. Wetness (= 1 for each wet face or 10 if flooded)
      varLabel(11+iColor) = 'wet_'//trim(geometry        &
        &                     %attribute                 &
        &                     %uni_name(sdr_seed_object) &
        &                     %val(iColor))
        
    end do    

    ! Set pointers appropriately.
    call sdr_append_protoVar( protoVar    = varLabel,   &
      &                       varSys      = varsys,     &
      &                       method_data = c_null_ptr, &
      &                       get_element = get_element )

    ! Create varMap
    call tem_create_varMap( varName = varsys%varname%val(       &
      &                               1:varsys%varname%nVals ), &
      &                     varSys  = varSys,                   &
      &                     varMap  = restart%varMap            )

    restart%header%headerPrefix = trim(main_debug%debugMesh)//trim(prefix_loc)&
      &                         //trim(charWave) //  '_restart_ProtoData'
    ! set the timestamp to be the current wave iWave stored in prefix
    restart%header%timestamp = trim(charWave)
    restart%header%binName = trim(main_debug%debugMesh) // trim(prefix_loc) &
      &                      // trim(charWave) // '_restart_ProtoData'      &
      &                      // EndianSuffix
    restart%write_file%nDofs = 1

    ! set the right solvertag
    restart%header%solverTag = tem_solverTag(header%general%solver)

    ! copy the global communicator to the restart one
    restart%comm = header%general%proc

    ! Dump global mesh info to header.lua.
    call dump_tem_global(tem_global)

    ! Write restart header information.
    call tem_restart_writeHeader( me     = restart,         &
      &                           tree   = dummy_treelmesh, &
      &                           varSys = varSys,          &
      &                           timing = timing           )

    ElemFileName = trim(tem_global%dirname) // 'elemlist' // EndianSuffix

    meshunit = newunit()

    inquire(iolength = rl) 0_long_k, 0_long_k

    open(unit = meshunit, file = trim(ElemFileName), action = 'write', &
      &  form = 'unformatted', access = 'direct', &
      &  recl = rl)

    restunit = newunit()

    inquire(iolength = rl) protodata

    open(unit = restunit, file = trim(restart%header%binName), &
      &  action = 'write', &
      &  form = 'unformatted', access = 'direct', &
      &  recl = rl)

    if (sdr_nodeProp_btest(node=proto%node, iNode=1, bit=isLeaf_bit)) then
      ! Only the root node exists as leaf.
      write(meshunit, rec=1) 0_long_k, 0_long_k
      write(restunit, rec=1) protoData_ofNode(1, proto, geometry, 0)
    else
      iLeaf = 0
      call write_childLeaves(meshunit = meshunit, restunit = restunit, &
        &                    iLeaf = iLeaf, level = 0, &
        &                    node_pos = 1, proto = proto, geometry = geometry)
    end if

    close(meshunit)
    close(restunit)

  end subroutine sdr_write_proto_as_restart
  ! ****************************************************************************


  ! ****************************************************************************
  !> Small helping routine to write leaves in order into a treelmesh formatted
  !! file.
  recursive subroutine write_childLeaves(meshunit, restunit, iLeaf, node_pos, &
    &                                    proto, geometry, level)
    integer, intent(in) :: meshunit
    integer, intent(in) :: restunit
    integer, intent(inout) :: iLeaf
    integer, intent(in) :: node_pos
    integer, intent(in) :: level
    type(sdr_prototree_type), intent(in) :: proto
    type(sdr_geometry_type), intent(in) :: geometry

    integer :: iChild
    integer :: childpos
    integer :: nextlevel

    nextlevel = level+1
    if( nextlevel > proto%nLevels ) return
    do iChild=0,7
      childpos = proto%node%linkpos(1)%val(node_pos) + iChild
      if (sdr_nodeProp_btest( node  = proto%node, &
        &                     iNode = childpos,   &
        &                     bit   = isLeaf_bit  ) ) then
        iLeaf = iLeaf + 1
        write(meshunit, rec=iLeaf) proto%node%treeID%val(childPos), 0_long_k
        write(restunit, rec=iLeaf) protoData_ofNode(childPos, proto, geometry, &
          &                                         level)
      else
        call write_childLeaves(meshunit = meshunit, restunit = restunit, &
          &                    iLeaf = iLeaf, &
          &                    node_pos = childpos, proto = proto, &
          &                    geometry = geometry, level = nextlevel)
      end if
    end do
  end subroutine write_childLeaves
  ! ****************************************************************************


  ! ****************************************************************************
  !> Small helping routine to get the variable data from a leaf.
  function protoData_ofNode(node_pos, proto, geometry, level) result(protodata)
    ! --------------------------------------------------------------------------!
    integer, intent(in) :: node_pos
    type(sdr_prototree_type), intent(in) :: proto
    type(sdr_geometry_type), intent(in) :: geometry
    integer, intent(in) :: level
    real(kind=rk) :: protodata(11+proto%node%nColors)
    ! --------------------------------------------------------------------------!

    integer :: iSide, iObj
    integer :: trumping, maxlevel
    integer :: obj_pos, attr_pos, attr_kind
    integer :: objkind_count(sdr_object_kinds_max)
    integer :: nObjects
    integer :: nColors
    integer :: iColor
    integer :: col_int, col_bit
    logical :: countAll
    integer :: intersected_first, intersected_last
    ! --------------------------------------------------------------------------!

    protodata = 0.0_rk
    objkind_count = 0
    nColors = proto%node%nColors

    if ( level==0 ) then
      nObjects = geometry%spatialObj%nVals
      intersected_first = 1
      intersected_last = nObjects
      countAll = .true.
    else
      intersected_first = proto%node%userObjPos%val(node_pos)%first
      intersected_last = proto%node%userObjPos%val(node_pos)%last
      nObjects = intersected_last - intersected_first + 1
      countAll = .false.
    end if

    ! - Refinement level
    protodata(1) = level

    ! - Wetness (= 1 for each wet face or 10 if flooded)
    ! Check flooding for first color (in first integer, bit 1):
    if ( sdr_nodeProp_btest(proto%node, node_pos, isFlooded_bit) ) then
      protodata(2) = 1.0_rk
    end if

    do iColor=1,proto%node%nColors
      col_int = (iColor-1) / proto%node%bytes_per_int + 1
      col_bit = mod(iColor-1, proto%node%bytes_per_int)*8 + 1
      if (btest(proto%node%PropertyBits%val(col_int,node_pos), col_bit)) then
        protodata(11+iColor) = 10.0_rk
      else
        do iSide=1,6
          if ( btest(proto%node%PropertyBits%val(col_int,node_pos), &
            &  col_bit+iSide) ) then
            protodata(11+iColor) = protodata(10+iColor) + 1.0_rk
          end if
        end do
      end if
    end do

    ! - Number of intersected objects
    protodata(3) = nObjects

    ! count number of each attibute kind
    maxlevel = 0
    trumping = geometry%attribute%dynArray%nVals + 1
    do iObj=intersected_first, intersected_last
      if (countAll) then
        obj_pos = iObj
      else
        obj_pos = proto%node%intersected_object%val(iObj) 
      end if
      attr_pos = geometry%spatialObj%val(obj_pos)%attribute_position
      attr_kind = geometry%attribute%dynArray%val(attr_pos)%kind
      objkind_count(attr_kind) = objkind_count(attr_kind) + 1
      if (attr_kind == sdr_Boundary_object) then
        trumping = min(trumping, attr_pos)
      end if
      maxlevel = max(maxlevel, geometry%attribute%dynArray%val(attr_pos)%level)
    end do

    ! - Intersected seeds
    protodata(4) = objkind_count(sdr_Seed_object)

    ! - Intersected refinements
    protodata(5) = objkind_count(sdr_Refinement_object)

    ! - Intersected fluidifyable elements
    protodata(6) = objkind_count(sdr_fluidifyable_object)

    ! - Intersected boundaries
    protodata(7) = objkind_count(sdr_Boundary_object)

    ! - Trumping boundary (trumping boundary condition or 0 if none)
    if (sdr_Boundary_object > 0) then
      protodata(8) = trumping
    end if

    ! - Target ref level (maximum level of intersected objects)
    protodata(9) = maxlevel

    ! - bc attribute position
    protodata(10) = proto%node%minBCID%val(node_pos)

    ! - hasBoundary bit is also set for non-flooded and intersected
    ! boundary bit so set hasBnd =1 only for flooded and hasboundary nodes
    if ( sdr_nodeProp_btest(proto%node, node_pos, hasBoundary_bit) .and. &
      & sdr_nodeProp_btest(proto%node, node_pos, isFlooded_bit) ) then
      protodata(11) = 1.0_rk
    end if

  end function protoData_ofNode
  ! ****************************************************************************


  ! ****************************************************************************
  !> Find the neighbor position in protoTree for iDir on the same level 
  !! or on any one above.
  function sdr_neighbor_in_proto(proto, coord, iDir, neighbor_level ) &
    & result( pos )
    !> preliminary tree
    type(sdr_protoTree_type), intent(in) :: proto
    integer, intent(in) :: coord(4)
    integer, intent(in) :: iDir
    integer, intent(out) :: neighbor_level
    integer :: pos
    ! --------------------------------------------------------------------------!
    integer :: myLevel, iLevel, offset(4), neighbor_coord(4)
    integer(kind=long_k) :: neighbor_tID

    pos = -1
    offset = 0
    neighbor_coord = 0
    myLevel = coord(4)
    offset(1:3) = qOffset(iDir,:)
    neighbor_coord = coord + offset
    neighbor_tID = tem_IdOfCoord( coord = neighbor_coord )
    ! Find the neighbor on the current level or on any one above.
    levelLoop: do iLevel = myLevel,0,-1
      pos = PositionOfVal(                              & 
        &     me = proto%node%treeID,                   &
        &     val = neighbor_tID,                       &
        &     lower = proto%levelNode_first( iLevel ),  &
        &     upper = proto%levelNode_last(  iLevel )   )
      if( pos > 0 ) then
        neighbor_level = iLevel
        return
      else
        ! Can NOT find neighbor, get its parent ID on one above level
        neighbor_tID = tem_ParentOf( neighbor_tID )
      end if
    end do levelLoop
  end function sdr_neighbor_in_proto
  ! ****************************************************************************

end module sdr_protoTree_module
