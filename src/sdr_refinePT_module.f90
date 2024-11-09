! Copyright (c) 2015-2016, 2019 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2015-2017, 2020 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2015 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
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
! ***************************************************************************
!> author: Kannan Masilamani
!! This module contains routine to refine protoTree until the minlevel
!! or level defined in the refinement object is reached
!!
module sdr_refinePT_module
  use env_module,            only: long_k, rk, globalMaxLevels, my_status_int
  use tem_param_module,      only: qQQQ, childCoord, qInvDir, qOffset
  use tem_dyn_array_module,  only: PositionOfVal, SortedPosOfVal, append
  use tem_grow_array_module, only: destroy, append, truncate, init, &
    &                              grw_longArray_type, grw_intArray_type
  use tem_logging_module,    only: logunit
  use tem_topology_module,   only: tem_directChildren, tem_FirstIdAtLevel,     &
    &                              tem_parentOf, tem_CoordOfId, tem_IdOfCoord
  use tem_timer_module,      only: tem_startTimer, tem_stopTimer 
  use tem_geometry_module,   only: tem_eligibleChildren
  use tem_tools_module,      only: tem_horizontalSpacer
  use tem_cube_module,       only: tem_cube_type
  use tem_sphere_module,     only: tem_sphere_type, tem_sphereCubeOverlap

  use sdr_spatialObj_module, only: sphere
  use sdr_prototree_module,  only: sdr_protoTree_type, levelValues_type,       &
    &                              sdr_neighbor_in_proto, sdr_node_neighbors
  use sdr_geometry_module,   only: sdr_geometry_type, is_intersecting
  use sdr_node_module,       only: intersectsBoundary_bit,  hasBoundary_bit,   &
    &                              isLeaf_bit, isTarget_bit, isFlooded_bit,    &
    &                              isFluidifyable_bit, append, truncate,       &
    &                              sdr_set_nodeProp_bit,                       &
    &                              sdr_clear_nodeProp_bit, sdr_nodeProp_btest, &
    &                              sdr_inHeritBnd_eligibleChildren,            &
    &                              sdr_append_childIntersectedObject,          &
    &                              sdr_intersectObjPos_type,                   &
    &                              grw_intersectObjPosArray_type
  use sdr_config_module,     only: sdr_confHead_type
  use sdr_timer_module,      only: timer_handle_refineLeaf,                    &
    &                              timer_handle_smoothLeaf,                    &
    &                              timer_handle_inHeritDR
  use sdr_attribute_module,  only: sdr_object_kinds_max,                       &
    &                              sdr_Fluidifyable_Object


  implicit none

  private

  public :: sdr_inHerit_distanceRefineObject
  public :: sdr_refine_leaf
  public :: sdr_smooth_leaf

contains
  
  ! ***************************************************************************
  !> This routines inherit distance refine sphere object from root node
  !! down to leaf node. 
  !! Only the object with level greater than node level are inHerited
  subroutine sdr_inHerit_distanceRefineObject(proto, geometry)
    !--------------------------------------------------------------------------!
    !> The proto tree description with all the data to refine further
    type(sdr_protoTree_type), intent(inout) :: proto
    !> type which contains all geometry object infos
    type(sdr_geometry_type), intent(in) :: geometry
    !--------------------------------------------------------------------------!
    integer :: iLevel, iParent, iChild
    integer :: firstParent, lastParent
    integer(kind=long_k) :: parent_ID_offset
    integer(kind=long_k) :: parentID
    integer :: parentProps
    !intersected distance objPos of parent node
    type(sdr_intersectObjPos_type) :: distObjPos
    type(levelValues_type) :: leVal
    logical :: testAll
    integer :: child_nodePos(8)
    integer :: nDistRefine_objs
    integer :: memLeft
    !--------------------------------------------------------------------------!
    call tem_startTimer( timerHandle = timer_handle_inHeritDR  ) 

    call tem_horizontalSpacer(funit=logunit(1))
    write(logUnit(5),*) 'Inherit distance refine objects ...'
    write(logUnit(5),*) 'Nr. of distance refine spatial objects: ', &
      &                 geometry%spatialObj%nVals - geometry%nUserObjs

    ! Set the length of the level cube length to the complete bounding cube
    ! length for level=0.
    ! Set some auxilary data describing the current level.
    leVal%dx = geometry%universe%extent
    leVal%ID_offset = tem_FirstIdAtLevel(0) ! first treeID on this level
    leVal%level = 0 ! level count

    testAll = .true.

    do iLevel = 1, proto%nLevels
      write(logunit(2),"(A,I0)") 'current Level: ', iLevel
      ! Position of first and last parent in the sorted treeID list.
      firstParent = proto%levelNode_first(iLevel-1)
      lastParent = proto%levelNode_last(iLevel-1)
      parent_ID_offset = leVal%ID_offset ! First treeID on parent level
      ! (save from previous iteration)

      ! Set some auxilary data describing the current level.
      leVal%dx = 0.5_rk * leVal%dx ! length of cube nodes on this level
      leVal%ID_offset = tem_FirstIdAtLevel(iLevel) ! first treeID on this level
      leVal%level = iLevel ! level count


      ! Iterate over all elements on the previous level.
      do iParent = firstParent, lastParent
        ! KM: Since no new leaf nodes are created, treeID array is already
        ! sorted so no need to look for sorted position

        ! get property of current node
        parentProps = ibits(proto%node%PropertyBits                    &
          &                      %val(proto%node%propLength, iParent), &
          &               pos = proto%node%lastbyte_pos, len=8         )
   
        ! Inherit only if parent is flooded and not a leaf
        if ( btest(parentProps, isFlooded_bit)  .and. &
          & .not. btest(parentProps, isLeaf_bit) ) then
   
          ! treeID of parent
          parentID = proto%node%treeID%val(iParent)

          if (testAll) then
            distObjPos%first = geometry%nUserObjs + 1
            distObjPos%last  = geometry%spatialObj%nVals
          else
            distObjPos%first = proto%node%distObjPos%val(iParent)%first
            distObjPos%last = proto%node%distObjPos%val(iParent)%last
          end if
          nDistRefine_objs = distObjPos%last - distObjPos%first + 1

          ! Children node position in protoTree
          child_nodePos(1) = proto%node%linkPos(1)%val(iParent)
          do iChild=2,8
            child_nodePos(iChild) = child_nodePos(1) + iChild - 1
          end do

          ! Inherit interested distance refine object to children
          ! and set distance_first and distance_last for each child
          call inHerit_intersectedObject(                                 &
            & proto                    = proto,                           &
            & geometry                 = geometry,                        &
            & parent                   = iParent,                         &
            & parentID                 = parentID,                        & 
            & testAll                  = testAll,                         &
            & intersected_object       = proto%node%intersected_distance, &
            & grwObjPos                = proto%node%distObjPos,           &
            & parent_ObjPos            = distObjPos,                      &
            & leVal                    = leVal,                           &
            & child_intersected_object = proto%child_intersected_object,  &
            & child_nodePos            = child_nodePos,                   &
            & memLeft                  = memLeft,                         &
            & isDistRefObj             = .true.                           )

          proto%node%memLeft_distObj = proto%node%memLeft_distObj + memLeft

        end if ! if not leaf
      end do !iNode
      testAll = .false.
      write(logUnit(5),"(A,I0)") 'Number of intersected distance objects:', &
        &                        proto%node%intersected_distance%nVals
    end do !iLevel

    write(logUnit(5),"(A,I0)") 'Number of intersected distance objects:', &
      &                  proto%node%intersected_distance%nVals

    write(logunit(10),*) 'Memory unused in parent intersected_distance: ', &
      &                 proto%node%memLeft_distObj
    call tem_stopTimer( timerHandle = timer_handle_inHeritDR )

  end subroutine sdr_inHerit_distanceRefineObject
  ! ***************************************************************************

  
  ! ***************************************************************************
  !> This routine extends the protoTree with max of minlevel or level of
  !! refinement object.
  !!
  !! If it is a leaf, check for intersected objects, and keep
  !! on refining accordingly down to the maximum requested level.
  !! When the desired level is reached (no refinement object,
  !! or desired level of intersected refinement object reached, check all
  !! (26) direct neigbors. If there is a neighbor intersected by a boundary
  !! check its refinement level, if it is higher then the current one,
  !! do another refinement step.
  subroutine sdr_refine_leaf(proto, geometry)
    !--------------------------------------------------------------------------!
    !> The proto tree description with all the data to refine further
    type(sdr_protoTree_type), intent(inout) :: proto
    !> type which contains all geometry object infos
    type(sdr_geometry_type), intent(in) :: geometry
    !--------------------------------------------------------------------------!
    integer(kind=long_k) :: parent_ID_offset
    integer(kind=long_k) :: parentID
    integer :: objkind_count(sdr_object_kinds_max)
    integer :: maxLevel
    integer :: parentProps
    integer :: iObject, iDist
    integer :: nObjects
    integer :: obj_pos, attr_pos, dist_pos
    integer :: attr_kind
    integer :: parent_pos
    integer :: iLevel, iParent
    integer :: firstParent, lastParent
    integer :: nNodesOld, nNodes
    logical :: testAll
    type(levelValues_type) :: leVal
    type(grw_longArray_type) :: grwTreeID
    !intersected user objects objPos of parent node
    type(sdr_intersectObjPos_type) :: userObjPos 
    !intersected distance objPos of parent node
    type(sdr_intersectObjPos_type) :: distObjPos
    integer :: child_nodePos(8)
    integer :: nDistRefine_objs
    integer :: memLeft
    logical :: noSub
    !integer :: VMHWM
    !--------------------------------------------------------------------------!
    call tem_startTimer( timerHandle = timer_handle_refineLeaf ) 

    call tem_horizontalSpacer(funit=logunit(1))
    write(logunit(1),*) 'Refining preliminary tree ...'

    ! Set the length of the level cube length to the complete bounding cube
    ! length for level=0.
    ! Set some auxilary data describing the current level.
    leVal%dx = geometry%universe%extent
    leVal%ID_offset = tem_FirstIdAtLevel(0) ! first treeID on this level
    leVal%level = 0 ! level count
 
    ! test all spatial objects in the root node
    testAll = .true.

    ! loop over all levels and refine flooded non-intersected leaf node to
    ! max of minlevel or level of refinement object
    levelLoop: do iLevel = 1, globalMaxlevels
      write(logunit(2),"(A,I0)") 'current level: ', iLevel

      ! Reset nNodes
      nNodes = 0

      ! Position of first and last parent in the sorted treeID list.
      firstParent = proto%levelNode_first(iLevel-1)
      lastParent = proto%levelNode_last(iLevel-1)
      parent_ID_offset = leVal%ID_offset ! First treeID on parent level
      ! (save from previous iteration)

      ! Initialize temporary growing array of treeID for current level
      call init(grwTreeID, length=lastParent-firstParent+1)

      ! Set some auxilary data describing the current level.
      leVal%dx = 0.5_rk * leVal%dx ! length of cube nodes on this level
      leVal%ID_offset = tem_FirstIdAtLevel(iLevel) ! first treeID on this level
      leVal%level = iLevel ! level count

      ! Current number of nodes in current level
      nNodesOld = proto%levelNode_last(iLevel)    &
        &       - proto%levelNode_first(iLevel) + 1

      ! Iterate over all elements on the previous level.
      parElemLoop: do iParent = firstParent, lastParent

        ! position of parent node in the sorted list
        parent_pos = proto%node%treeID%sorted(iParent)
        objkind_count = 0

        ! get property of parent node 
        parentProps = ibits(proto%node%PropertyBits                       &
          &                      %val(proto%node%propLength, parent_Pos), &
          &               pos = proto%node%lastbyte_pos, len=8            )

        ! if parent is leaf, flooded and non-intersected by boundary
        ! than check for maxLevel of the parent for further refinement
        if ( btest(parentProps, isFlooded_bit) .and. &
          &  btest(parentProps, isLeaf_bit)  ) then 

          ! Only need to refine further if this node is not intersected by a
          ! boundary and not a target node or node with subelement resolution.
          ! As they always are fully resolved during the
          ! protoTree generation before.
          if (proto%node%subelement_resolution>0) then
            noSub = proto%node%subLevel%val(parent_pos) < 0
          else
            noSub = .true.
          end if

          if ( .not. btest(parentProps, intersectsBoundary_bit) .and. &
            & (noSub .or. btest(parentProps, isTarget_bit)) ) then

            ! Initialize maxLevel
            maxLevel = proto%node%minLevel%val(parent_pos)
        
            ! treeID of parent
            parentID = proto%node%treeID%val(parent_pos)

            ! test sphere object for distance refinement and get max refine
            ! level
            if (proto%node%distanceRefine) then
              distObjPos%first = proto%node%distObjPos%val(parent_pos)%first
              distObjPos%last = proto%node%distObjPos%val(parent_pos)%last
              nDistRefine_objs = distObjPos%last - distObjPos%first + 1
              do iDist = distObjPos%first, distObjPos%last
                dist_pos = proto%node%intersected_distance%val(iDist)
                attr_pos = geometry%spatialObj%val(dist_pos)%attribute_position
                maxLevel = max(maxLevel, &
                  &            geometry%attribute%dynArray%val(attr_pos)%level)
              end do
            else
              nDistRefine_objs = 0
            end if
             
            ! Step1: Get max refine level among all intersected objects of the
            !        parent.
            if (testAll) then
              userObjPos%first = 1
              userObjPos%last = geometry%nUserObjs
            else
              userObjPos%first = proto%node%userObjPos%val(parent_pos)%first
              userObjPos%last = proto%node%userObjPos%val(parent_pos)%last
            end if
            nObjects = userObjPos%last - userObjPos%first + 1
            do iObject = userObjPos%first, userObjPos%last
              if (testAll) then
                obj_pos = iObject
              else
                obj_pos = proto%node%intersected_object%val(iObject) 
              end if
              attr_pos = geometry%spatialObj%val(obj_pos)%attribute_position
              attr_kind = geometry%attribute%dynArray%val(attr_pos)%kind
              objkind_count(attr_kind) = objkind_count(attr_kind) + 1
              maxLevel = max(maxLevel, &
                &            geometry%attribute%dynArray%val(attr_pos)%level)
            end do

            ! Maximum level is already reached, if there are fluidifyable objects
            ! then set isFluidifyable_bit
            if ( maxLevel < leVal%level .and. &
              & objkind_count(sdr_Fluidifyable_object) > 0 ) then
              call sdr_set_nodeProp_bit( node  = proto%node,        &
                &                        iNode = parent_pos,        &
                &                        bit   = isFluidifyable_bit )
            end if

            ! Step2: if smoothBounds are true, get level of neighbor 
            ! which is intersected with boundary
            if (geometry%smoothbounds .and. &
              & btest(parentProps, hasBoundary_bit)) then
              call check_bndLevel( proto            = proto,            &
                &                  parent           = parent_pos,       &
                &                  parent_ID_offset = parent_ID_offset, &
                &                  leVal            = leVal,            &
                &                  maxLevel         = maxLevel          )
            end if

            ! Step3: If current node is not at its finest level or a neighbor
            ! node with intersected_boundary is on a finer level, refine the
            ! current node:
            ! Append 8 children in the protoTree
            notreachedmax: if ( maxLevel >= leVal%level ) then
              ! Create 8 childrens and inherit property bits from parent
              call create_children( proto         = proto,         &
                &                   parent        = parent_pos,    &
                &                   child_nodePos = child_nodePos, &
                &                   grwTreeID     = grwTreeID      )
              
              ! Inherit intersected objects from parent to child
              call inHerit_intersectedObject(                                  &
                & proto                    = proto,                            &
                & geometry                 = geometry,                         &
                & parent                   = parent_pos,                       &
                & parentID                 = parentID,                         &
                & testAll                  = testAll,                          &
                & intersected_object       = proto%node%intersected_object,    &
                & grwObjPos                = proto%node%userObjPos,            &
                & parent_objPos            = userObjPos,                       &
                & leVal                    = leVal,                            &
                & child_intersected_object = proto%child_intersected_object,   &
                & child_nodePos            = child_nodePos,                    &
                & memLeft                  = memLeft,                          &
                & isDistRefObj             = .false.                           )

              proto%node%memLeft_userObj = proto%node%memLeft_userObj + memLeft  

              ! Inherit intersected distance refine objects from parent to child
              if (proto%node%distanceRefine) then
                call inHerit_intersectedObject(                                &
                  & proto                    = proto,                          &
                  & geometry                 = geometry,                       &
                  & parent                   = parent_pos,                     &
                  & parentID                 = parentID,                       &
                  & testAll                  = .false.,                        &
                  & intersected_object       = proto%node%intersected_distance,&
                  & grwObjPos                = proto%node%distObjPos,          &
                  & parent_ObjPos            = distObjPos,                     &
                  & leVal                    = leVal,                          &
                  & child_intersected_object = proto%child_intersected_object, &
                  & child_nodePos            = child_nodePos,                  &
                  & memLeft                  = memLeft,                        &
                  & isDistRefObj             = .true.                          )

                proto%node%memLeft_distObj = proto%node%memLeft_distObj + memLeft  
              end if

              ! update number of nodes created in this level
              nNodes = nNodes + 8
            end if notreachedmax
          end if ! not intersected boundary
        end if ! flooded, leaf bit
      end do parElemLoop
  
      ! Copy treeID of new leaf nodes created in this level to dynamic array
      ! of treeIDs
      call append(proto%node%treeID, grwTreeID%val(:nNodes))
      ! destroy temporary growing array of treeID
      call destroy(grwTreeID)

      ! Set levelNode last for current level with number of nodes 
      ! created in current level
      proto%levelNode_last(iLevel) = proto%levelNode_last(iLevel) + nNodes

      ! update levelNode first and last for next levels   
      proto%levelNode_first(iLevel+1:) = proto%levelNode_first(iLevel+1:) &
        &                              + nNodes
      proto%levelNode_last(iLevel+1:) = proto%levelNode_last(iLevel+1:) &
        &                             + nNodes

      ! The number of elements on the current level
      write(logUnit(5),"(A,I0)") '  *   new nodes on this level: ', nNodes
      write(logUnit(5),"(A,I0)") '  * Total nodes on this level: ', nNodes + nNodesOld
      write(logUnit(5),"(A)")    ''

      !VMHWM =  my_status_int('VmHWM:')
      !write(logUnit(3),*) 'Memory usuage:', VMHWM

      ! Leave the loop, if no new nodes were created in this iteration.
      if (proto%levelNode_first(iLevel) > proto%levelNode_last(iLevel)) exit

      ! Update the number of levels in the tree, as the current level contains
      ! new elements and level higher than previous nLevels
      if (proto%nLevels < iLevel) proto%nLevels = iLevel

      ! after root node deactivate testAll
      testAll = .false.

    end do levelLoop

    ! destroy temporary intersected object
    call destroy(proto%child_intersected_object)

    ! Clean up memory and free those entries that are not actually needed.
    call truncate(proto%node)

    call tem_stopTimer( timerHandle = timer_handle_refineLeaf )

    write(logunit(10),*) 'Memory unused in parent intersected_object: ', &
      &                 proto%node%memLeft_userObj

    write(logunit(10),*) 'Memory unused in parent intersected_distance: ', &
      &                 proto%node%memLeft_distObj

  end subroutine sdr_refine_leaf
  ! ***************************************************************************


  ! ***************************************************************************
  !> This routine inherit the intersected boundary objects from parent to
  !! childrens
  !!
  !! Update the intersected_first and last for children 
  subroutine inHerit_intersectedObject(proto, geometry, parent, parentID,      &
    &                                  testAll, intersected_object,            &
    &                                  grwObjPos, parent_objPos,               &
    &                                  leVal, child_intersected_object,        &
    &                                  child_nodePos, memLeft, isDistRefObj    )
    !--------------------------------------------------------------------------!
    !> The proto tree description with all the data to refine further
    type(sdr_protoTree_type), intent(inout) :: proto
    !> type which contains all geometry object infos
    type(sdr_geometry_type), intent(in) :: geometry
    !> Position of parent node in protoTree
    integer, intent(in) :: parent
    !> TreeID of parent
    integer(kind=long_k), intent(in) :: parentID
    !> To test all intersected objects
    logical, intent(in) :: testAll
    !> Growing array of intersected objects. 
    !! Could be user defined  or distance refine spatial objects
    type(grw_intArray_type), intent(inout) :: intersected_object
    !> First and last position of intersected object of all nodes in
    !! intersected_object list
    type(grw_intersectObjPosArray_type), intent(inout) :: grwObjPos
    !> Position of first ans last intersected object of parent node 
    !! in intersected_object list
    type(sdr_intersectObjPos_type), intent(in) :: parent_objPos
    !> contains information on current level on which children are created
    type(levelValues_type), intent(in) :: leVal
    !> Temporary array of intersected objects for 8 children
    type(grw_intArray_type), intent(inout) :: child_intersected_object
    !> 8 children node position in protoTree
    integer, intent(in) :: child_nodePos(8)
    !> memory of parent intersected object unused by children
    integer, intent(out) :: memLeft
    !> Is this distance refine objects
    logical, intent(in) :: isDistRefObj
    !--------------------------------------------------------------------------!
    integer(kind=long_k) :: firstchild_ID
    integer :: node_coord(4)
    integer :: firstchild_coord(4)
    integer :: iChild, iObject, obj_pos
    integer :: child_nObjects(8)
    integer :: objLevel
    integer :: attr_pos
    integer :: minLevel
    integer :: prim_pos
    integer :: prim_kind
    type(sdr_intersectObjPos_type) :: child_objPos(8)

    type(tem_cube_type) :: node_cube
    type(tem_sphere_type) :: hol_sphere ! hollow sphere
    !--------------------------------------------------------------------------!

    ! Get first treeID of parent
    firstchild_ID = parentID * 8_long_k + 1_long_k

    ! Get the length of children cube. all children cubes have same length.
    node_cube%extent = leVal%dx
    node_cube%halfwidth = 0.5_rk * leVal%dx

    ! Coordinate of childs can be computed with an offset from
    ! first child using childCoord parameter
    firstchild_coord = tem_coordOfId(treeID = firstchild_ID, &
      &                              offset = leVal%ID_offset)

    ! Initialize counter for child_intersected_first and last     
    child_objPos(:)%first = 1
    child_objPos(:)%last = 0
    child_intersected_object%nVals = 0

    ! Step5: create children and inherit the intersected object information
    childLoop: do iChild = 1, 8
      ! minimum level of node
      minLevel = proto%node%minLevel%val(parent)
    
      ! Define the children cube to intersect with.
      node_coord(1:3) = firstchild_coord(1:3) + childCoord(iChild,:) 
      node_cube%origin = leVal%dx * node_coord(1:3) + geometry%universe%origin
      node_cube%endPnt = leVal%dx * (node_coord(1:3)+1) &
        &                + geometry%universe%origin
      node_cube%center = node_cube%origin + node_cube%halfwidth

      ! the first intersected object of this child has to be one after
      ! all the previous childs intersected objects
      child_objPos(iChild:)%first = child_intersected_object%nVals + 1
      ! Loop over the objects which intersected with parent to 
      ! determine objects intersected with child
      do iObject = parent_objPos%first, parent_objPos%last
        if (testAll) then
          obj_pos = iObject
        else
          obj_pos = intersected_object%val(iObject)
        end if
        attr_pos = geometry%spatialObj%val(obj_pos)%attribute_position
        objLevel = geometry%attribute%dynArray%val(attr_pos)%level

        ! Test for intersection of geometric object and node cube.
        if (is_intersecting(node_cube, geometry, obj_pos)) then
          ! Inherit objects only if objLevel >= current level and minLevel
          if (objLevel > max(minLevel, leVal%level))  then
            ! For sphere objects check if node is inside sphere
            ! of intersected with sphere surface.
            ! This can be determined by checking intersection with
            ! hollow sphere
            prim_pos = geometry%spatialObj%val(obj_pos)%primitive_position
            prim_kind = geometry%spatialObj%val(obj_pos)%geometry_primitive
            if (prim_kind==sphere .and. isDistRefObj) then
              hol_sphere = geometry%sphere%val(prim_pos)
              hol_sphere%only_surface = .true.
              ! if hollow sphere intersect with node then append
              ! this object to intersected object list
              if (tem_sphereCubeOverlap(hol_sphere, node_cube)) then
                call append( child_intersected_object, obj_pos )
              else
                ! node is inside sphere
                !
                ! maximum level amoung intersected objects 
                minLevel = max(minLevel, objLevel)
              end if
            else  
              ! If geometry is not a sphere then add the object to the
              ! intersected_object list of the child.
              call append( child_intersected_object, obj_pos )

            end if !if object is sphere
          end if ! objLevel>mylevel
        end if !intersected
      end do !iObject
      
      ! Set minLevel to which this node to be refined
      proto%node%minLevel%val(child_nodePos(iChild)) = minLevel

      ! Set child intersected object first and last in child_intersected_object
      child_objPos(iChild:)%last = proto%child_intersected_object%nVals
      child_nObjects(iChild) = child_objPos(iChild)%last    &
        &                    - child_objPos(iChild)%first + 1 
    end do childLoop

    ! Now copy temporary growing array of intersected objects of child
    ! to actual array of intersected_objects
    call sdr_append_childIntersectedObject(                              &
      &             geometry                 = geometry,                 &
      &             node                     = proto%node,               &
      &             parent                   = parent,                   &
      &             testAll                  = testAll,                  &
      &             intersected_object       = intersected_object,       &
      &             grwObjPos                = grwObjPos,                &
      &             child_nodePos            = child_nodePos,            &
      &             child_intersected_object = child_intersected_object, &
      &             child_objPos             = child_objPos,             &
      &             memLeft                  = memLeft                   )

  end subroutine inHerit_intersectedObject
  ! ***************************************************************************


  ! ***************************************************************************
  !> This routine checks if neighbor node with intersected boundary is level 
  !! higher than current node level.
  !!
  !! If neighbor node is intersected boundary bit but no a leaf or target or
  !! node with qVal than refine current node.
  subroutine check_bndLevel(proto, parent, parent_ID_offset, leVal, maxLevel)
    !--------------------------------------------------------------------------!
    !> preliminary tree on which childern are created
    type(sdr_protoTree_type), intent(inout) :: proto
    !> Position of parent node on the dynamic array of node%treeID and node_data
    !! in preliminary tree
    integer, intent(in) :: parent
    !> first treeID of the parent
    integer(kind=long_k), intent(in) :: parent_ID_offset
    !> contains information on current level on which children are created
    type(levelValues_type), intent(in) :: leVal
    !> Maximum level to refine current node
    integer, intent(inout) :: maxLevel
    !--------------------------------------------------------------------------!
    integer :: iDir
    integer :: coord(4)
    integer :: neighbor_level
    integer :: neighbor_pos
    logical :: refine_for_neighbor
    !--------------------------------------------------------------------------!
    ! Get coordinate of current element
    coord = tem_CoordOfId( treeID = proto%node%treeID%val(parent), &
      &                    offset = parent_ID_offset               )

    ! loop over all 26 directions
    do iDir = 1, qQQQ
      ! get position of neighbor in the protoTree which might me
      ! from different level
      neighbor_pos = sdr_neighbor_in_proto( proto, coord, iDir,  &
        &                                   neighbor_level )
        
      if ( sdr_nodeProp_btest(node  = proto%node,           &
        &                     iNode = neighbor_pos,         &
        &                     bit   = intersectsBoundary_bit) ) then
        ! Neighbor is intersected by a boundary, now we check if we need to
        ! refine further down.

        ! If neighbor_pos is not a leaf or target node, refine the fluid
        ! element to bnd_level if smoothbounds are active.
        refine_for_neighbor &
          &  = .not. ( sdr_nodeProp_btest(node   = proto%node,         &
          &                               iNode  = neighbor_pos,       &
          &                               bit    = isLeaf_bit    )     &
          &            .or. sdr_nodeProp_btest(node  = proto%node,     &
          &                                    iNode = neighbor_pos,   &
          &                                    bit   = isTarget_bit  ) )

        if (refine_for_neighbor) then
          ! Only refine, if smoothbounds are requested, or the neighbor is
          ! the actual boundary. This might happen, if the boundary is
          ! covering a complete element on a coarser level, then it is
          ! refined to and the flooded neighbor is not refined down to that
          ! level yet.
          ! This affects all parent levels and therefore might influence a
          ! large part of the mesh, which might not be desired, thus it can
          ! be deactivated by the smoothbounds setting.
          maxLevel = max(maxLevel, leVal%level)
        end if
      end if  
    end do

  end subroutine check_bndLevel
  ! ***************************************************************************


  ! ***************************************************************************
  !> This routine smoothens fluid domain with maximum level jumps of 1.
  !!
  !! This is done by running over all leaf, flooded, non-intersected node
  !! at each level and check for neighbor in all 26 directions.
  !! If neighbor exist in protoTree and is not a leaf then get eligible
  !! children of neighbor in the inverse direction of iDir and check 
  !! if any of eligible children is not a leaf then refine myself.
  !!
  !! Algorithm:
  !! Iterate over minLevel, maxLevel-2
  !! - Iterate over all nodes in iLevel
  !!   + if iNode is leaf, flooded and non intersected boundary
  !!     * Iterate over 26 directions
  !!       ++ if neighbor in iDir exist in protoTree and not a leaf
  !!          ** get eligible children of neighbor in inverse iDir
  !!            -- if any of eligible children in protoTree is not a leaf
  !!               +++ refine iNode
  subroutine sdr_smooth_leaf(proto, header, maxLevel)
    !--------------------------------------------------------------------------!
    !> preliminary tree on which childern are created
    type(sdr_protoTree_type), intent(inout) :: proto
    !> some global information on solver name and version
    type(sdr_confHead_type), intent(inout) :: header
    !> Maximum level in the fluid domain
    integer, intent(in) :: maxLevel
    !--------------------------------------------------------------------------!
    integer :: iLevel, iWave, iParent, iDir, iChild
    integer :: parent_pos, neighbor_pos
    integer :: parentProps 
    integer :: nNodes_total, nNodes_level, nNodesOld
    integer :: firstParent, lastParent
    integer :: coord(4)
    integer :: offset(4), neighbor_coord(4)
    integer(kind=long_k) :: neighbor_tID
    integer(kind=long_k) :: parent_ID_offset
    integer, allocatable :: eligible_childs(:)
    integer :: eligible_childPos
    logical :: isChildLeaf(8)
    logical :: isChildTarget(8)
    type(grw_longArray_type) :: grwTreeID
    integer :: child_nodePos(8)
    logical :: noSub
    logical :: neigh_noLeaf, neigh_noTarget, neigh_noSub
    !--------------------------------------------------------------------------!
    call tem_startTimer( timerHandle = timer_handle_smoothLeaf ) 

    call tem_horizontalSpacer(funit=logunit(1))
    write(logunit(1),*) 'Smooth levels ...'
    iWave = 0
    ! Leave loop when no new nodes are created 
    waves: do 
      ! Wave count
      iWave = iWave + 1
      write(logunit(2),"(A,I0)") ' current wave ', iWave
      ! number of nodes in protoTree
      nNodes_total = proto%node%nNodes

      ! Start from 2nd coarsest level to 2nd finest level.
      ! In every wave, we could reduce loop to -2 level since
      ! all neccessary nodes are already created in that level.
      ! In this loop, we create new nodes of the iLevel
      ! so we run over nodes in iLevel-1 to create children in iLevel
      do iLevel = header%minLevel+1, maxLevel - iWave
        write(logunit(2),"(A,I0)") 'current level ', iLevel
  
        ! Position of first and last treeID in current level in sorted treeID list
        firstParent = proto%levelNode_first(iLevel-1)
        lastParent = proto%levelNode_last(iLevel-1)
        ! first treeID on parent level
        parent_ID_offset = tem_FirstIdAtLevel(iLevel-1)
 
        ! Initialize temporary growing array of treeID for current level
        call init(grwTreeID, length=lastParent-firstParent+1)

        ! Count number of new nodes created in this level
        nNodes_level = 0

        ! Current number of nodes in current level
        nNodesOld = proto%levelNode_last(iLevel)    &
          &       - proto%levelNode_first(iLevel) + 1

        ! Iterate over all nodes in previous level to create nodes 
        ! in current level
        do iParent = firstParent, lastParent
  
          ! position of parent node in the sorted list
          parent_pos = proto%node%treeID%sorted(iParent)
  
          ! get property of parent node 
          parentProps = ibits(proto%node%PropertyBits                     &
            &                    %val(proto%node%propLength, parent_pos), &
            &               pos = proto%node%lastbyte_pos, len=8        )
  
          ! if parent node is leaf, flooded  and non-intersected
          ! by boundary than check for neighbor level
          if ( btest(parentProps, isFlooded_bit) .and. &
            &  btest(parentProps, isLeaf_bit)  ) then 

            ! Only need to refine further if this node is not intersected by a
            ! boundary and not a target node or node with subelement resolution.
            ! As they always are fully resolved during the
            ! protoTree generation before.
            if (proto%node%subelement_resolution>0) then
              noSub = proto%node%subLevel%val(parent_pos) < 0
            else
              noSub = .true.
            end if

            if ( .not. btest(parentProps, intersectsBoundary_bit) .and. &
              & (noSub .or. btest(parentProps, isTarget_bit)) ) then

              ! Get coordinate of current element
              coord = tem_CoordOfId(                                  &
                &         treeID = proto%node%treeID%val(parent_pos), &
                &         offset = parent_ID_offset                   )
  
              offset = 0
              neighbor_coord = 0
              ! loop over all 26 directions and get position of neighbor in 
              ! protoTree
              dirLoop: do iDir = 1, qQQQ
                offset(1:3) = qOffset(iDir,:)
                neighbor_coord = coord + offset
                neighbor_tID = tem_IdOfCoord( coord = neighbor_coord )
                neighbor_pos = PositionOfVal(                               &
                  &              me    = proto%node%treeID,                 & 
                  &              val   = neighbor_tID,                      &
                  &              lower = proto%levelNode_first( iLevel-1 ), &
                  &              upper = proto%levelNode_last(  iLevel-1 )  )
  
                ! neighbor tID exist in protoTree on the iLevel
                if (neighbor_pos > 0) then
                  ! create children when neighbor is not a leaf and
                  ! neither target bit or subelement
                  neigh_noLeaf = .not. sdr_nodeProp_btest(node  = proto%node,   &
                    &                                     iNode = neighbor_pos, &
                    &                                     bit   = isLeaf_bit    ) 
                  
                  neigh_noTarget = .not. sdr_nodeProp_btest(node  = proto%node,  &
                    &                                       iNode = neighbor_pos,&
                    &                                       bit   = isTarget_bit ) 

                  if (proto%node%subelement_resolution>0) then
                    neigh_noSub = proto%node%subLevel%val(neighbor_pos) < 0
                  else
                    neigh_noSub = .true.
                  end if

                  ! if neighbor is virtual node
                  if (neigh_noLeaf .and. (neigh_noSub .or. neigh_noTarget)) then
                    ! treeIDs of children
                    ! childID(1:8) = tem_directChildren(neighbor_tID)
                    ! get eligible children in inverse of iDir
                    call tem_eligibleChildren(eligible_childs, qInvDir(iDir))

                    ! set default all child to leaf
                    isChildLeaf = .true.
                    isChildTarget = .true.
                    do iChild = 1, size(eligible_childs)
                      ! position of eligible childrens in protoTree
                      ! since they are continous in memory we could find 
                      ! it with an offset given by eligible_childs
                      ! and linkPos(1) of parent
                      eligible_childPos = proto%node%linkPos(1)%val(neighbor_pos) &
                        &               + eligible_childs(iChild) - 1
                      isChildLeaf((eligible_childs(iChild))) =            &
                        &  sdr_nodeProp_btest( node  = proto%node,        &
                        &                      iNode = eligible_childPos, &
                        &                      bit   = isLeaf_bit         ) 

                      isChildTarget((eligible_childs(iChild))) =          &
                        &  sdr_nodeProp_btest( node  = proto%node,        &
                        &                      iNode = eligible_childPos, &
                        &                      bit   = isTarget_bit       ) 
                    end do

                    deallocate(eligible_childs)
                    ! if any child is not a leaf or target refine my self
                    if ( .not. all(isChildLeaf .or. isChildTarget) ) then
                      !write(*,*) 'isChildLeaf ', isChildLeaf
                      ! For this childrens, we inherit only property from
                      ! parent. Interected objects are not inherited since
                      ! parent it already refined at maximum level of its 
                      ! refinement object.
                      call create_children( proto         = proto,         &
                        &                   parent        = parent_pos,    & 
                        &                   child_nodePos = child_nodePos, &
                        &                   grwTreeID     = grwTreeID      )
                      
                      nNodes_level = nNodes_level + 8
                      ! exit dirLoop and go to next node in level
                      exit dirLoop
                    end if ! atleast one child is not leaf
                  end if !neighbor is not a leaf
                end if !neighbor exist in protoTree
              end do dirLoop
            end if ! not intersected boundary
          end if ! flooded leaf node
        end do ! iNode
  
        ! Copy treeID of new leaf nodes created in this level to dynamic array
        ! of treeIDs
        call append(proto%node%treeID, grwTreeID%val(:nNodes_level))
        ! destroy temporary growing array of treeID
        call destroy(grwTreeID)

        ! Set levelNode last for next level with number of nodes 
        ! created in next level
        proto%levelNode_last(iLevel) = proto%levelNode_last(iLevel) &
           &                           + nNodes_level
      
        ! update levelNode first and last for next levels   
        proto%levelNode_first(iLevel+1:) = proto%levelNode_first(iLevel+1:) &
          &                              + nNodes_level
        proto%levelNode_last(iLevel+1:) = proto%levelNode_last(iLevel+1:) &
          &                             + nNodes_level

        ! The number of elements on the current level
        write(logUnit(3),"(A,I0)") '  *   new nodes on this level: ', nNodes_level
        write(logUnit(3),"(A,I0)") '  * Total nodes on this level: ', &
          &                        nNodes_level + nNodesOld
        write(logunit(3),*) ''

      end do ! ilevel
      
      ! exit loop when no new nodes are created
      if (nNodes_total == proto%node%nNodes) exit waves

    end do waves
    write(logUnit(1),"(A,I0)") ' done in nWaves ', iWave

    call tem_stopTimer( timerHandle = timer_handle_smoothLeaf )
  end subroutine sdr_smooth_leaf
  ! ***************************************************************************


  ! ***************************************************************************
  !> This routine append 8 children to protoTree and inherit property bits from 
  !! parent. leaf bit is removed from parent.
  subroutine create_children(proto, parent, child_nodePos, grwTreeID)
    !--------------------------------------------------------------------------!
    !> preliminary tree on which childern are created
    type(sdr_protoTree_type), intent(inout) :: proto
    !> Position of parent node on the dynamic array of node%treeID and node_data
    !! in preliminary tree
    integer, intent(in) :: parent
    !> 8 children node position in protoTree
    integer, intent(out) :: child_nodePos(8)
    !> Temporary growing array of TreeID contains new leaf nodes in current level
    type(grw_longArray_type), intent(inout) :: grwTreeID
    !--------------------------------------------------------------------------!
    integer :: iChild
    integer(kind=long_k) :: treeID(1:8)
    integer :: child_sublevel
    integer :: propbits(proto%node%propLength)
    logical :: child_hasBnd(8)
    integer :: minLevel
    !--------------------------------------------------------------------------!
    ! If parent has hasBoundary_bit then inherit this property to its 
    ! eligible children
    child_hasBnd = sdr_inHeritBnd_eligibleChildren(proto%node, parent)

    ! Get treeID of eight children
    treeID(1:8) = tem_directChildren(proto%node%treeID%val(parent))
    
    ! Inherit the property bits to children 
    propBits = proto%node%PropertyBits%val(:, parent)

    ! remove leaf bit from parent node
    call sdr_clear_nodeProp_bit( node  = proto%node, &
      &                          iNode = parent,     &
      &                          bit   = isLeaf_bit  )

    ! set position of first child at linkPos(1) of parent node
    proto%node%linkPos(1)%val(parent) = proto%node%nNodes + 1

    ! Increase the counter for proto%nLeafNodes and proto%nFloodedLeaves
    ! Since we create 8 new children in place of 1 parent we just
    ! increase this counter by 7
    proto%nLeafNodes = proto%nLeafNodes + 7
    proto%nFloodedLeaves= proto%nFloodedLeaves + 7

    if (proto%node%subelement_resolution > 0) then
      child_sublevel = proto%node%sublevel%val(parent) - 1
    else
      ! If there is no subelement resolution the sublevels are not available.
      ! However, we store a negative value here in child_sublevel to avoid
      ! further checks for subelement_resolution further down.
      child_sublevel = -1
    end if

    minLevel = proto%node%minLevel%val(parent)
    ! create children and inherit the intersected object information
    childLoop: do iChild = 1, 8
      ! append this child to tree and inherit minLevel of parent to child
      ! since for distance refine objects if node is inside this 
      ! distance than refine this node to minLevel
      call append(me           = proto%node,                      &
        &         treeID       = treeID(iChild),                  &
        &         grwTreeID    = grwTreeID,                       &
        &         PropertyBits = propBits,                        &
        &         sublevel     = child_sublevel,                  &
        &         minLevel     = minLevel,                        &
        &         pos          = child_nodePos(iChild)            )

      ! if this child hasBnd, set the hasBoundary_bit 
      ! and copy the directions which has boundary from parent.
      ! if parent has boundary in certain direction, its child should have 
      ! the save since hasBoundary is set by the leaf of intersected 
      ! boundary node
      if (child_hasBnd(iChild)) then
        call sdr_set_nodeProp_bit(node  = proto%node,            &
          &                       iNode = child_nodePos(iChild), &
          &                       bit   = hasBoundary_bit        )
        proto%node%hasBndBits%val(child_nodePos(iChild)) =    &
          &                   proto%node%hasBndBits%val(parent)
      end if
    end do childLoop

  end subroutine create_children
  ! ***************************************************************************

end module sdr_refinePT_module
! ***************************************************************************
