! Copyright (c) 2012-2015, 2017, 2019 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2012, 2014 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2013-2014, 2023 Harald Klimach <harald.klimach@dlr.de>
! Copyright (c) 2014 Matthias Johannink
! Copyright (c) 2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2017 Raphael Haupt <Raphael.Haupt@student.uni-siegen.de>
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
!! This module contains routines for boundary identification
!! and qVal computation
!!
!!@todo HK: This module is utterly confusing, I do not really understand what
!!          all this code is actually doing, and it is hard to parse it.
!!          The module needs to be redesigned and cleaned up!
module sdr_boundary_module
  use env_module,            only: long_k, rk, newunit, isLittleEndian, &
    &                              globalMaxLevels
  use treelmesh_module,      only: treelmesh_type
  use tem_param_module,      only: qQQQ, qOffset, qInvDir
  use tem_mergesort_module,  only: mrgrnk
  use tem_topology_module,   only: tem_directChildren, tem_firstIDatLevel, &
    &                              tem_CoordOfId, tem_IdOfCoord,           &
    &                              tem_parentOf, tem_LevelOf
  use tem_geometry_module,   only: tem_eligibleChildren, tem_BaryOfCoord, &
    &                              tem_CoordOfReal, tem_BaryOfID,         &
    &                              tem_ElemSize, tem_vrtxCoordOfId
  use tem_line_module,       only: tem_line_type, intersect_RayTriangle, &
    &                              fraction_PointLine
  use tem_plane_module,      only: tem_plane_type
  use tem_cube_module,       only: tem_cube_type, tem_convertTreeIDtoCube
  use tem_point_module,      only: tem_point_type, tem_pointCubeOverlap
  use tem_triangle_module,   only: tem_triangle_close_face, &
    &                              tem_triangle_normal_proximity
  use tem_dyn_array_module,  only: PositionOfVal
  use tem_grow_array_module, only: grw_realArray_type, init, append, destroy
  use tem_aux_module,        only: tem_abort
  use tem_logging_module,    only: logunit
  use tem_debug_module,      only: dbgUnit

  use sdr_prototree_module,  only: sdr_protoTree_type, levelValues_type, &
    &                              sdr_neighbor_in_proto
  use sdr_geometry_module,   only: sdr_geometry_type
  use sdr_periodic_module,   only: sdr_periodicPlane_type
  use sdr_node_module,       only: sdr_node_type, isFlooded_bit, isLeaf_bit, &
    &                              isTarget_bit, IntersectsBoundary_bit,     &
    &                              isFluidifyable_bit, sdr_nodeProp_btest
  use sdr_spatialObj_module, only: periodicPlane, triangle
  use sdr_attribute_module,  only: sdr_attrList_type, sdr_Boundary_object

  implicit none

  private

  public :: sdr_identify_boundary
  public :: needCalcQValByBCID
  public :: needFldDglByBCID
  public :: sdr_qValByNode
  public :: sdr_find_periodic_neighbor

  !> default qVal when no intersection
  real(kind=rk), parameter, public :: sdr_qVal_no_intersection = -1._rk


contains


  ! ************************************************************************** !
  !> This routine checks for boundary neighbors and level of the boundary
  !! node
  !!
  !! Note, this can not easily be used for target nodes with subelement
  !! resolution, as it assumes q-Values if the node is intersected by a
  !! boundary.
  subroutine sdr_identify_boundary(node_pos, treeID, coord, leVal, proto, &
    &                              geometry, BC_ID, qVal, bc_normal,      &
    &                              meshUniverse                           )
    ! ---------------------------------------------------------------------- !
    !> Position of leaf in the preliminary tree
    !@todo: the position of what node?
    integer, intent(in) :: node_pos
    !> treeID of parent node
    integer(kind=long_k), intent(in) :: treeID
    !> Coordinate of treeID
    integer, intent(in) :: coord(4)
    !> level value of parent node
    type(levelValues_type), intent(in) :: leVal
    !> preliminary tree
    type(sdr_protoTree_type), intent(in) :: proto
    !> contains all geometrical objects
    type(sdr_geometry_type), intent(in) :: geometry
    !> Boundary ID for all 26 neighbor directions
    integer(kind=long_k), intent(out) :: BC_ID(qQQQ)
    !> distance from boundary for all 26 neighbor directions
    real(kind=rk), intent(out) :: qVal( qQQQ )
    !> Wall normal pointing from the surface to the barycenter of the
    !! boundary element
    real(kind=rk), intent(out) :: bc_normal(3)
    !> contains bounding cube information
    type(treelmesh_type), intent(in) :: meshUniverse
    ! ---------------------------------------------------------------------- !
    integer :: iDir, neighbor_pos
    integer :: neighbor_level
    logical :: unKnownBnd(qQQQ)
    integer :: startDir
    integer :: unKnownNeighPos(qQQQ)
    real(kind=rk) :: elembary(3)
    real(kind=rk) :: neigh_normal(3,qQQQ)
    real(kind=rk) :: neigh_distsq(qQQQ)
    integer(kind=long_k) :: minBCID
    logical :: check_unKnownBnd
    logical :: check_neigh_normals
    integer :: iDir_tmp
    integer :: nBCs
    integer :: nNeigh_normals
    integer(kind=long_k) :: neigh_bcid
    type(tem_cube_type) :: neighcube
    type(tem_point_type) :: neighpoint
    ! --------------------------------------------------------------------------!
    !write(dbgUnit(1),*) 'Entering identify boundary'
    !To detect hanging node in non-direct boundary direction with no
    !boundary attached to it
    unKnownBnd = .false.
    unKnownNeighPos = -1
    ! No need to check for unknownBnd if there is level jump
    ! i.e. bndLevel > level
    check_unKnownBnd = .true.
    ! compute minBCID of this node and use it for unKnownBnd
    minBCID = huge(minBCID)

    qVal = sdr_qVal_no_intersection
    bc_normal = 0.0_rk
    ! Get barycenter of current element
    elemBary = tem_BaryOfID( meshUniverse, treeID )
    !write(dbgUnit(1),*) 'myID: ', treeID
    !write(dbgUnit(1),*) 'myBary: ', elemBary

    ! First check if current node is intersecting a boundary
    ! then qVal must be active, as nodes with subelement resolution are
    ! treated seperately.
    ! So check for intersection of
    ! link with the geometry intersected by this node.
    ! If there is no intersection then check its neighbor
    if ( sdr_nodeProp_btest(node  = proto%node,           &
      &                     iNode = node_pos,             &
      &                     bit   = intersectsBoundary_bit) ) then

      if ( .not. sdr_nodeProp_btest(node  = proto%node,       &
        &                           iNode = node_pos,         &
        &                           bit   = isFluidifyable_bit) ) then
        linkloop: do iDir = 1,qQQQ
          call getBCID_and_calcQval( proto        = proto,            &
            &                        geometry     = geometry,         &
            &                        elemBary     = elemBary,         &
            &                        iDir         = iDir,             &
            &                        bndnode_pos  = node_pos,         &
            &                        level        = leVal%level,      &
            &                        leVal        = leVal,            &
            &                        meshUniverse = meshUniverse,     &
            &                        bc_id        = BC_ID(iDir),      &
            &                        minBCID      = minBCID,          &
            &                        unKnownBnd   = unKnownBnd(iDir), &
            &                        qVal         = qVal(iDir)        )
        end do linkloop
      end if

      bc_normal = getNormal( geometry, elemBary, proto%node, node_pos, &
        &                    int(minBCID)                              )

    end if

    check_neigh_normals = (dot_product(bc_normal, bc_normal) < tiny(bc_normal(1)))

    if (check_neigh_normals) then
      nNeigh_normals = 0
    end if

    ! For nodes which have a boundary neighbor
    do iDir = 1, qQQQ

      ! Is the current element with boundaries intersected in iDir?
      noq: if ( qVal(iDir) == sdr_qVal_no_intersection ) then

        ! get position of neighbor in the protoTree which might me
        ! from different level
        neighbor_pos = sdr_neighbor_in_proto( proto, coord, iDir, &
          &                                   neighbor_level      )

        unKnownNeighPos(iDir) = neighbor_pos

        bndBit: if ( sdr_nodeProp_btest(node  = proto%node,           &
          &                             iNode = neighbor_pos,         &
          &                             bit   = intersectsBoundary_bit) ) then

          ! Neighbor is intersected by a boundary look up the actual boundary
          ! condition if the neighbor needs qvalue computations or is not
          ! flooded.
          !
          ! Only consider neighbors that are not flooded by any color for
          ! boundary conditions or that require qvalues, if the neighbor is
          ! flooded, this is an internal boundary that probably only separates
          ! two different colors.
          ! For those we do not want to create boundary conditions unless we
          ! have to compute qValues and there is geometry between the current
          ! element and the flooded, intersected element.
          neigh_bcid = int(proto%node%minBCID%val(neighbor_pos), kind=long_k)
          if ( .not. sdr_nodeProp_btest(node  = proto%node,                  &
            &                           iNode = neighbor_pos,                &
            &                           bit   = isFlooded_bit )              &
            & .or. needCalcQValByBCID( geometry%attribute, int(neigh_bcid) ) &
            & ) then
            call getBCID_and_calcQval( proto        = proto,            &
              &                        geometry     = geometry,         &
              &                        elemBary     = elemBary,         &
              &                        iDir         = iDir,             &
              &                        bndnode_pos  = neighbor_pos,     &
              &                        level        = neighbor_level,   &
              &                        leVal        = leVal,            &
              &                        meshUniverse = meshUniverse,     &
              &                        bc_id        = BC_ID(iDir),      &
              &                        minBCID      = minBCID,          &
              &                        unKnownBnd   = unKnownBnd(iDir), &
              &                        qVal         = qVal(iDir)        )
          end if

          if (check_neigh_normals) then
            neigh_normal(:, nNeigh_normals+1) = getNormal(                   &
              &   geometry, elemBary, proto%node, neighbor_pos, int(minBCID) )
            neighpoint%coord = elemBary - neigh_normal(:, nNeigh_normals+1)
            call tem_convertTreeIDtoCube(cube = neighcube,                &
              &                          tree = meshUniverse,             &
              &                          treeID = proto%node              &
              &                                        %treeID            &
              &                                        %val(neighbor_pos) )
            if (tem_pointCubeOverlap(neighpoint, neighcube)) then
              nNeigh_normals = nNeigh_normals + 1
              neigh_distsq(nNeigh_normals) &
                &  = dot_product(neigh_normal(:,nNeigh_normals), &
                &                neigh_normal(:,nNeigh_normals)  )
            end if
          end if

        else bndBit
          if ( .not. sdr_nodeProp_btest(node  = proto%node,   &
            &                           iNode = neighbor_pos, &
            &                           bit   = isFlooded_bit ) ) then
            ! This non-flooded element does not contain a boundary condition.
            ! Usually this should not happen, however due to additional
            ! flooding of qValue elements it might be, that a neighbor is in
            ! the non-flooded domain.
            ! This might also happen if there is an element with subelement
            ! resolution adjacent to a non-flooded domain.
            ! Special treatment is needed for boundaries in this direction.
            unKnownBnd(iDir) = .true.
            unKnownNeighPos(iDir) = neighbor_pos
          end if
        end if bndBit

      end if noq

      !write(dbgUnit(5),*) 'iDir = ', iDir, '; bcID = ', bc_ID(iDir),  &
      !  &                 '; Q-Value = ',qVal(iDir)
    end do ! iDir

    if (check_neigh_normals .and. (nNeigh_normals > 0)) then
      bc_normal = neigh_normal(:, minloc(neigh_distsq(:nNeigh_normals), 1))
    end if

    ! set min bcid for element with no boundary id
    ! set qVal = 0.5
    nBCs = geometry%attribute%kindpos(sdr_Boundary_object)%nVals
    checkUnknowns: if (check_unKnownBnd .and. any(unKnownBnd)) then
      do iDir=1,qQQQ
        unknownDir: if (unKnownBnd(iDir)) then
          bc_id(iDir) = minBCID
          write(logUnit(3),*) iDir, 'minBCID of unKnownBnd: ', minBCID
!          bc_id(iDir) = minval(bc_id, bc_id >0)
          ! need calc qVal?
          if( needCalcQValByBCID( geometry%attribute, int(bc_id(iDir)) )) &
            & qVal(iDir) = 0.5_rk

          if (bc_id(iDir) == 0 .or. bc_id(iDir) > nBCs) then
            write(logunit(0),*) 'ERROR: Unable to treat the undefined node'
            write(logunit(0),*) '       No valid boundary IDs are found'
            write(logunit(0),*) 'startDir    :', startDir
            write(logunit(0),*) 'iDir        :', iDir
            write(logunit(0),*) 'minBCID     :', minBCID
            write(logunit(0),*) 'qVal        :', qVal(iDir)
            write(logunit(0),*) 'treeID      :', treeID
            write(logunit(0),*) 'neighbor_pos:', unknownNeighPos(iDir)
            write(logunit(0),*) 'neigh_bary  :',                          &
              &             tem_baryOfID( meshUniverse,                   &
              &                           proto%node%treeID               &
              &                                %val(unknownNeighPos(iDir)))
            write(logunit(0),*) 'propbit     :',                 &
              &            proto%node%PropertyBits               &
              &                      %val(:,unknownNeighPos(iDir))
            write(logUnit(1),*) 'Other directions:'
            do iDir_tmp=1,qQQQ
              if (iDir_tmp /= iDir) then
                write(logUnit(1),*) 'iDir: ', iDir_tmp
                write(logUnit(1),*) 'Bits:',                           &
                  &            proto%node%PropertyBits                 &
                  &                      %val(:,unknownNeighPos(iDir_tmp))
                write(logUnit(1),*) 'neigh_bary  :',                      &
                  &             tem_baryOfID( meshUniverse,               &
                  &                           proto%node%treeID           &
                  &                                %val(unknownNeighPos(  &
                  &                                       iDir_tmp)))
              end if
            end do
            call tem_abort()
          end if
        end if unknownDir
      end do
    end if checkUnknowns
    !write(dbgUnit(5),*) &
    !  &     '-----------  ROUTINE OUT: sdr_identify_boundary   -----------'
  end subroutine sdr_identify_boundary
  ! ****************************************************************************


  ! ****************************************************************************
  !> This routine gets minBCID of the given node position in the protoTree.
  !! If the minBcid is periodic then it bcID is set to treeID of fluid node
  !! on the opposite side of periodic plane.
  !! It also computes the qVal if calc_dist = true. If qVal = -1 then
  !! there is no intersection and if qVal > 1 then the geometry is intersected
  !! after the link distance.
  subroutine getBCID_and_calcQval(proto, geometry, elemBary, iDir,         &
    &                             bndnode_pos, level, leVal, meshUniverse, &
    &                             bc_id, minBCID, qVal, unKnownBnd )
    ! --------------------------------------------------------------------------!
    !> preliminary tree
    type(sdr_protoTree_type), intent(in) :: proto
    !> contains all geometrical objects
    type(sdr_geometry_type), intent(in) :: geometry
    !> current element barycenter
    real(kind=rk), intent(in) :: elemBary(3)
    !> Current boundary neighbor direction
    integer, intent(in) :: iDir
    !> position of node treeID in the proto%node list
    integer, intent(in) :: bndnode_pos
    !> level of the node
    integer, intent(in) :: level
    !> level value of parent node
    type(levelValues_type), intent(in) :: leVal
    !> contains bounding cube information
    type(treelmesh_type), intent(in) :: meshUniverse
    !> Boundary ID for direction iDir
    integer(kind=long_k), intent(out) :: BC_ID
    !> minimum boundary id of current node before truncation
    integer(kind=long_k), intent(inout) :: minBCID
    !> distance from boundary for direction iDir
    real(kind=rk), intent(out) :: qVal
    !> Is true if a neighbor with no property is encountered
    logical, intent(inout) :: unKnownBnd
    ! --------------------------------------------------------------------------!
    !write(dbgUnit(1),*) 'iDir ', iDir
    bc_id = 0_long_k
    qVal = sdr_qVal_no_intersection

    ! Get BC ID for this link
    BC_ID = int(proto%node%minBCID%val(bndnode_pos), kind=long_k)

    !!write(dbgUnit(5),*) 'min bc_id ', bc_id

    ! If BC_ID == -1 then it is periodic
    ! because periodic boundary id is set to -1
    if (bc_id == -1_long_k) then
      call sdr_find_periodic_neighbor( elemBary       = elemBary,    &
        &                              iDir           = iDir,        &
        &                              bc_id          = BC_ID,       &
        &                              qVal           = qVal,        &
        &                              unKnownBnd     = unKnownBnd,  &
        &                              neighbor_pos   = bndnode_pos, &
        &                              neighbor_level = level,       &
        &                              leVal          = leVal,       &
        &                              proto          = proto,       &
        &                              geometry       = geometry,    &
        &                              meshUniverse   = meshUniverse )
    end if ! periodic boundary

    ! set minBCID before truncation
    ! needed to set unknown boundary
    ! KM: bc_id >0 - no periodic boundary
    if (BC_ID>0)  minBCID = min(minBCID, BC_ID)

    ! JQ: calculate q-values
    ! in some cases, might qVal > 1 or qVal == -1.0,
    ! take special treatment in sdr_truncate_qVal
    ! KM: periodic boundary might have encountered node with qVal on
    ! periodic domain and it might have already computed qVal.
    ! So do compute qVal only if there is no_intersection
    if ( needCalcQValByBCID( geometry%attribute, int(BC_ID) ) &
      & .and. qVal == sdr_qVal_no_intersection ) then
      call sdr_qValByNode( proto, geometry, leVal%dx, iDir, &
        & elemBary, bndnode_pos, qVal )
      call sdr_truncate_qVal( proto    = proto,      &
        &                     qVal     = qVal,       &
        &                     BCID     = bc_id,      &
        &                     neighPos = bndnode_Pos )
    endif

    !write(dbgUnit(5),*) 'final bc_id ', bc_id
    !write(dbgUnit(5),*) 'final qVal ', qVal

  end subroutine getBCID_and_calcQval
  ! ****************************************************************************


  ! ****************************************************************************
  !> This routine find the treeID on the opposite side neighbor of the
  !! periodic plane for current leaf node
  subroutine sdr_find_periodic_neighbor( elemBary, iDir, bc_id, qVal,          &
    &                                    unKnownBnd, neighbor_pos,             &
    &                                    neighbor_level, leVal, proto,         &
    &                                    geometry, meshUniverse )
    ! --------------------------------------------------------------------------!
    !> current element barycenter
    real(kind=rk), intent(in) :: elemBary(3)
    !> Current boundary neighbor direction
    integer, intent(in) :: iDir
    !> treeiD of opposite neighbor with periodic plane
    integer(kind=long_k), intent(inout) :: bc_ID

    !> distance from boundary for direction iDir
    !!
    !! @todo HK: maybe turn this optional, why is it relevant for periodic?
    real(kind=rk), intent(out) :: qVal

    !> Set to true if a neighbor with no property is encountered
    logical, intent(inout) :: unKnownBnd

    !> preliminary tree
    type(sdr_protoTree_type), intent(in) :: proto
    !> contains all geometrical objects
    type(sdr_geometry_type), intent(in) :: geometry
    !> level of the periodic boundary neighbor node
    integer, intent(in) :: neighbor_level
    !> position of neighbor treeID in the proto%node list
    integer, intent(in) :: neighbor_pos
    !> level value of parent node
    type(levelValues_type), intent(in) :: leVal
    !> contains bounding cube information
    type(treelmesh_type), intent(in) :: meshUniverse
    ! --------------------------------------------------------------------------!
    integer :: opp_periplane_pos, periPlane_pos
    type( sdr_periodicPlane_type ) :: plane_curr, plane_opp
    real(kind=rk) :: bary_opp(3)
    real(kind=rk) :: coordReal_opp_per(3)
    !> position of current periodic plane in periodic plane list
    integer :: iObj
    integer :: coord_opp(4)
    real(kind=rk) :: bary_opp_per(3)
    integer(kind=long_k) :: treeID_opp
    integer :: treeID_opp_pos
    integer :: nBCs
    integer :: minbcid, leafLevel
    logical :: bc_defined
    integer :: intersected_first, intersected_last
    integer :: minLevel_loc
    ! --------------------------------------------------------------------------!
    minLevel_loc = min(leVal%level, neighbor_level)
    !write(dbgUnit(1),*) 'myLevel: ', leVal%level, ' neighborLevel: ', neighbor_level
    qVal = sdr_qVal_no_intersection

    periPlane_pos = 0
    opp_periplane_pos = 0

    !starting coordinate for projecting is current element barycenter
    coordReal_opp_per = elemBary

    intersected_first = proto%node%userObjPos%val(neighbor_pos)%first
    intersected_last = proto%node%userObjPos%val(neighbor_pos)%last
    ! get the periodoc plane position in the array of periodic planes
    do iObj = intersected_first, intersected_last
      if(geometry%spatialObj&
        &        %val( proto%node%intersected_object &
        &                        %val(iObj) )%geometry_primitive &
        & .eq. periodicPlane ) then

        periPlane_pos = geometry%spatialObj &
          &                     %val( proto%node%intersected_object &
          &                                     %val(iObj) )%primitive_position
        plane_curr = geometry%periPlane%val( periPlane_pos )

        opp_periPlane_pos = geometry%spatialObj%val( &
          &                 geometry%periPlane%val( &
          &                 periPlane_pos )%oppPlane_pos)%primitive_position

        plane_opp = geometry%periPlane%val( opp_periPlane_pos )

        ! Project current element on current periodic plane
        ! and translate to opposite plane
        coordReal_opp_per = projectVecOnPlane( &
          &                (coordReal_opp_per - plane_curr%plane%origin), &
          &                                     plane_curr%plane ) + &
          &                plane_opp%plane%origin
      endif
    enddo  !periodic loop

    ! barycenter of opposite element in the requested direction iDir from
    ! the opposite periodic plane must be a fluid element
    ! we are looking for. Check for existence of this fluid element
    ! in the protoTree later before assigning to bcid

    ! barycenter of element in opposite periodic plane
    bary_opp_per = tem_BaryOfID( meshUniverse, &
      &               tem_IdOfCoord( &
      &                  tem_CoordOfReal( meshUniverse, &
      &                      coordReal_opp_per, leVal%level ) ) )

    ! coordinate of element in opposite periodic plane in fluid level
    coord_opp = tem_coordOfID( &
      &           tem_IdOfCoord( &
      &              tem_CoordOfReal( meshUniverse, &
      &                  bary_opp_per, leVal%level ) ) )

    ! Get the coordinate of element in the direction of boundary search.
    ! This element is suppose to be fluid element
    coord_opp(1:3) = coord_opp(1:3) + qOffset(iDir,:)

    ! treeID of opposite suspected fluid element
    treeID_opp = tem_IdOfCoord( coord_opp )
    !write(dbgUnit(5),*) 'treeID_opp ', treeID_opp

    ! Barycenter of opposite suspected fluid element
    bary_opp = tem_BaryOfID( meshUniverse, treeID_opp )
    !write(dbgUnit(5),*) "bary_opp  : ",bary_opp

    ! Get the position of treeID of bary_opp on the current level since the
    ! solver requires neighbor treeID in the same level .
    ! Also the treeID must be a leaf because its property are checked later
    treeID_opp_pos = getTreeIDPosOfCoord( coordReal  = bary_opp,        &
      &                                   mesh       = meshUniverse,    &
      &                                   minLevel   = minLevel_loc,    &
      &                                   maxLevel   = globalMaxLevels, &
      &                                   leafLevel  = leafLevel,       &
      &                                   leafTreeID = treeID_opp,      &
      &                                   proto      = proto            )

    !write(dbgUnit(5),*) 'treeID_opp_pos ', treeID_opp_pos
    !write(dbgUnit(5),*) 'leaf level, treeID: ', leafLevel, treeID_opp

    ! if projected suspected fluid element exist in protoTree
    ! check if min bcid is again periodic.
    ! This happens if periodic plane cuts the element in the barycenter
    ! and both elements in the fine level has periodic boundary property.
    if(treeID_opp_pos > 0) then
      ! get min BCID from intersected_objects of current code
      bc_id = proto%node%minBCID%val(treeID_opp_pos)
      !write(dbgUnit(1),*) '1st assign bc_id ', bc_id
      ! if bc_id == -1, its again periodic.
      ! happens when periodic plane passes through barycenter of the parent
      ! node
      if (bc_id == -1) then
        !write(dbgUnit(1),*) 'Its again periodic'
        ! do last three steps again
        ! ie. compute new bary_opp by offsetting in normal*dx
        ! and get the TreeID of bary_opp and treeID position in protoTree.
        ! Normal of current plane points outwards the fluid domain .
        bary_opp = bary_opp + plane_curr%plane%unitNormal * leVal%dx
        !write(dbgUnit(5),*) "bary_opp  : ",bary_opp

        ! get the treeID of bary_opp
        treeID_opp_pos = getTreeIDPosOfCoord( coordReal  = bary_opp,        &
          &                                   mesh       = meshUniverse,    &
          &                                   minLevel   = minLevel_loc,    &
          &                                   maxLevel   = globalMaxLevels, &
          &                                   leafTreeID = treeID_opp,      &
          &                                   leafLevel  = leafLevel,       &
          &                                   proto      = proto            )
      endif
    endif

    ! Projected fluid element must exist in protoTree since
    ! protoTree contains final fluid elements.
    ! If does not exist then something else is wrong
    bc_defined = .false.
    checkTID: if (treeID_opp_pos > 0 .and. leafLevel > 0) then


      ! if flooded and min bc_id == 0 then its fluid node
      ! without any intersecting boundary then return treeID of
      ! opposite side of periodic boundary on the level of fluid element
      if (sdr_nodeProp_btest( node  = proto%node,     &
        &                     iNode = treeID_opp_pos, &
        &                     bit   = isFlooded_bit   ) ) then
        bc_id = -tem_IdOfCoord( tem_CoordOfReal( meshUniverse,          &
          &                                      bary_opp, leVal%level) )
        bc_defined = .true.
        !write(dbgUnit(1),*) 'Is flooded'
        !write(dbgUnit(1),*) 'bc_id in fluid_level: ', bc_id
      else if ( sdr_nodeProp_btest(node  = proto%node,        &
        &                     iNode = treeID_opp_pos,         &
        &                     bit   = intersectsBoundary_bit) &
        &                                                     ) then
        ! if interesected and exist valid bc_id then return that bc_id
        ! in special case for qVal active for that bc_id, compute
        ! qVal and do special treatment in sdr_truncate_qVal.
        !write(dbgUnit(1),*) 'Is intersected boundary'

        ! get min BCID from intersected_objects of current code
        minbcid = proto%node%minBCID%val(treeID_opp_pos)
        !write(dbgUnit(1),*) 'minBCID: ', minbcid

        nBCs = geometry%attribute%kindpos(sdr_boundary_object)%nVals
        if (minbcid > 0 .and. minbcid <= nBCs) then
          bc_id = minbcid
          bc_defined = .true.
          ! if qVal active for this bc_id
          if( needCalcQValByBCID(geometry%attribute, minbcid) ) then
            call sdr_qValByNode( proto, geometry, leVal%dx, iDir, &
              &                  bary_opp_per, treeID_opp_pos, qVal )

            !write(dbgUnit(5),*) 'qVal ', qVal
            ! if qVal > 1 or no intersection and
            ! current node is also flooded then set bc_id = -treeID_opp
            ! else keep valid bc_id and set qVal = 1.0.
            call sdr_truncate_qVal( proto           = proto,          &
              &                     qVal            = qVal,           &
              &                     BCID            = bc_id,          &
              &                     neighPos        = treeID_opp_Pos, &
              &                     treeID_periodic = treeID_opp      )
            !!write(dbgUnit(5),*) '1.1 qVal ', qVal
          end if
          !write(dbgUnit(5),*) 'boundary bc_id', bc_id
          ! return valid bc_id if qVal < 1 else return negative treeID
          ! of periodic neighbor
          return
        end if
      else
        ! Periodic neighbor node is not flooded and not intersected i.e
        ! undefined hanging node, so treat it with minBCID from all direction
        ! of current node in  sdr_identify_boundary
        write(logUnit(3),*) 'WARNING: Periodic neighbor is undefined node.'
        write(logUnit(3),*) '         Setting minBCID from other directions!'
        write(logUnit(5),*) 'iDir: ', iDir, 'ElemBary: ', elemBary
        write(logUnit(5),*) 'bary_opp', bary_opp
        unKnownBnd = .true.
        return

      end if

    else checkTID

      ! projected fluid treeID not found in protoTree.
      ! It must be a treeID inside the boundary which has no property and
      ! also not refined to the same level as boundary so does not exist in
      ! protoTree
      bc_id = - treeID_opp
      unKnownBnd = .true.
      write(logunit(3),*) 'WARNING: Projected fluid element not found in &
        &protoTree'
      !write(logUnit(0),*) 'TreeID_opp_pos ', treeId_opp_pos, ' leaflvl:', leafLevel
      write(logUnit(5),*) 'iDir: ', iDir, 'ElemBary: ', elemBary
      write(logUnit(5),*) 'bary_opp', bary_opp
      return

    end if checkTID

    !write(dbgUnit(1),*) 'bcID: ', bc_id
    !flush(dbgUnit(1))
    if (.not. bc_defined) then
      write(logUnit(0),*) 'TreeID_opp_pos ', treeId_opp_pos, ' leaflvl:', leafLevel
      write(logUnit(0),*) 'property ', proto%node%propertyBits%val(:,treeID_opp_pos)
      ! If the boundary is not properly defined, the
      ! direction of the plane is wrong or
      ! the element has undefined property i.e not flooded nor
      ! intersected bit due to hollow geometries.
      write(logunit(0),*) 'ERROR: in periodic projection. Neighbor element of'
      write(logunit(0),*) ' the opposide periodic plane has undefined property.'
      write(logUnit(0),*) bary_opp
      write(logUnit(0),*) 'iDir: ', iDir, 'ElemBary: ', elemBary
      write(logunit(0),*) 'HINT 1: Check the normal direction of the plane' &
        &                 // ' ID', opp_periPlane_pos
      write(logunit(0),*) '  is wrong.'
      write(logunit(0),*) 'Solution: Swap the periodic plane vectors to' &
        &                 // ' change the normal.'
      write(logunit(0),*) '          Normal of both periodic must point' &
        &                 // ' opposite'
      write(logunit(0),*) '          fluid domain'
      write(logunit(0),*) 'HINT 2: If above solution does not work then '
      write(logunit(0),*) '        shift the plane position by dx/2'
      write(logunit(0),*) 'HINT 3: projected element might be inside the hollow'
      write(logunit(0),*) '        geometry. deactivate hollow geometry'
      write(logunit(0),*) 'HINT 4: If above solution does not help then' &
        &                 // ' periodic'
      write(logunit(0),*) '        boundary interface does not match.'
      write(logunit(0),*) '        check debug output for any mismatch on the'
      write(logunit(0),*) '        periodic boundary planes refinement'
      call tem_abort()
    end if

  end subroutine sdr_find_periodic_neighbor
  ! ****************************************************************************


  ! ****************************************************************************
  !> This function returns the position of treeID of given coordReal in the
  !! the given mesh
  !! Start from minLevel which is the level of neighbor and find the treeID
  !! which is a leaf in protoTree
  function getTreeIDPosOfCoord( coordReal, mesh, minLevel, maxLevel, &
    &                           leafLevel, leafTreeID, proto ) result(pos)
    ! --------------------------------------------------------------------------!
    real(kind=rk), intent(in) :: coordReal(3)
    !>Mesh contain geometry universe (bounding cube) info
    type( treelmesh_type ), intent(in) :: mesh
    ! treeID of given coordReal in bounding cube in leaf level
    integer(kind=long_k), intent(inout) :: leafTreeID
    !> minlevel
    integer, intent(in) :: minLevel
    !> maxlevel
    integer, intent(in) :: maxLevel
    !> level in which leaf node is found
    integer, intent(out) :: leafLevel
    !> preliminary tree
    type(sdr_protoTree_type), intent(in) :: proto
    !> position of treeID in proto tree
    integer :: pos
    ! --------------------------------------------------------------------------!
    integer :: iLevel
    ! --------------------------------------------------------------------------!
    leafLevel = 0
    pos = 0

    do iLevel = minLevel, maxLevel
      leafTreeID = tem_IdOfCoord( tem_CoordOfReal( mesh, coordReal, iLevel) )

      pos = PositionOfVal(                                    &
        &              me    = proto%node%treeID,             &
        &              val   = leafTreeID,                    &
        &              lower = proto%levelNode_first(iLevel), &
        &              upper = proto%levelNode_last(iLevel)   )

      ! found leaf return the current treeID position
      if (sdr_nodeProp_btest(node  = proto%node,    &
        &                    iNode = pos,           &
        &                    bit   = isLeaf_bit) ) then
        leafLevel = iLevel
        return
      end if
    end do

  end function getTreeIDPosOfCoord
  ! ****************************************************************************


  ! *****************************************************************************
  !> This function project given vector on an given plane
  !!
  !!Example: projection of vector a onto a vector u is given as
  !! \f$ proj_u a = \frac{a \cdot u}{|u|^2} \cdot u \f$
  ! *****************************************************************************
  function projectVecOnPlane( vecU, plane ) result(res)
  ! ---------------------------------------------------------------------------!
  !> vector to project
  real(kind=rk), intent(in) :: vecU(3)
  !> plane on which vecU will be projected
  type( tem_plane_type ), intent(in) :: plane
  !> output projected value
  real(kind=rk) :: res(3)
  ! ---------------------------------------------------------------------------!

  ! project vecU on plane%vecA and plane%vecB and add them up
  ! to find the projection of vecU on the plane
  res = ( dot_product( vecU, plane%vec(:,1)) /       &
    &     dot_product( plane%vec(:,1), plane%vec(:,1) ) ) * plane%vec(:,1) &
    &     + ( dot_product( vecU, plane%vec(:,2)) /   &
    &     dot_product( plane%vec(:,2), plane%vec(:,2) ) ) * plane%vec(:,2)

  end function projectVecOnPlane
  ! *****************************************************************************


  ! *****************************************************************************
  !> This routine computes the minimum distance of a given link and all the
  !! geometries in a given node:\n
  !! the link is given by a vector and a origin point.\n
  !! the node is given by the node position in the protoTree.\n
  !! If there is no intersection, qVal returns -1.0
  subroutine sdr_qValByNode(proto, geometry, dx, iDir, origin, &
    &                           node_pos, qVal)
    ! --------------------------------------------------------------------------!
    !> preliminary tree
    type(sdr_protoTree_type), intent(in) :: proto
    !> contains all geometrical objects
    type(sdr_geometry_type), intent(in) :: geometry
    !> dx of current level
    real(kind=rk), intent(in) :: dx
    !> Direction
    integer, intent(in) :: iDir
    !> position of node to find the geometries
    integer :: node_pos
    !> current element barycenter
    real(kind=rk), intent(in) :: origin(3)
    !> distance from boundary for all 26 neighbor directions
    real(kind=rk), intent(out)  :: qVal
    ! --------------------------------------------------------------------------!
    integer :: iSpatialObj
    integer :: obj_pos, geom_prim, atb_pos, prim_pos
    logical :: intersected
    type( tem_line_type ) :: line
    type( tem_point_type ) :: intersect_p
    real(kind=rk) :: qVal_t
    integer :: intersected_first, intersected_last
    ! --------------------------------------------------------------------------!
    qVal_t = huge( qVal_t )
    ! create line by origin, iDir and dx
    line%origin = origin
    line%vec = qOffset(iDir, :) * dx

    ! JQ: loop over spatialObjects, check whether need to calc dist
    ! call routine that calc distance, given vector and triangles
    ! if intersected, return the smallest distance
    ! if not intersected, return 0.0
    intersected_first = proto%node%userObjPos%val(node_pos)%first
    intersected_last = proto%node%userObjPos%val(node_pos)%last
    do iSpatialObj = intersected_first, intersected_last
      ! Object position
      obj_pos = proto%node%intersected_object%val(iSpatialObj)
      ! Geometry primitive
      geom_prim = geometry%spatialObj%val(obj_pos)%geometry_primitive
      ! geometry primitive position
      prim_pos = geometry%spatialObj%val(obj_pos)%primitive_position
      ! Attribute position
      atb_pos = geometry%spatialObj%val(obj_pos)%attribute_position

      ! check whether to calculate qVal
      if( geometry%attribute%dynArray%val(atb_pos)%calc_dist ) then
        select case( geom_prim )
          case( triangle )

          ! Compute whether the link (as described by line) is intersecting
          ! this triangle:
          intersected = intersect_RayTriangle(                          &
            &             line = line,                                  &
            &             triangle = geometry%triangle%val( prim_pos ), &
            &             intersect_p = intersect_p                     )
          if (intersected) then
            ! calculate fraction as qVal
            ! compare it with previous one, choose smaller one
            qVal_t = min( qVal_t, fraction_PointLine(intersect_p, line) )
          endif

          !@todo: add other geometry primitive here
          case default
        end select ! geom_prim
      endif ! whether to calculate qVal

    enddo ! spatial objects loop

    ! if ever intersected, take it as output; otherwise, return -1.0
    if ( qVal_t /= huge(qVal_t) ) then
      qVal = qVal_t
    else
      qVal = sdr_qVal_no_intersection
    endif

  end subroutine sdr_qValByNode
  ! ****************************************************************************

  ! ****************************************************************************!
  !> This routine checks if a boundary need calc qVal for a given BCID
  !! It is used in identify_boundary routine
  function needCalcQValByBCID( attribute, bcid ) result( calc_qVal )
    type(sdr_attrList_type), intent(in) :: attribute
    integer, intent(in) :: bcid
    integer :: iAtt
    logical :: calc_qVal

    if ((bcid <= attribute%kindpos(sdr_Boundary_object)%nVals) &
      & .and. (bcid >= 1)) then
      iAtt = attribute%kindpos(sdr_Boundary_object)%val(bcid)
      calc_qVal = attribute%dynArray%val(iAtt)%calc_dist
    else
      calc_qVal = .false.
    end if
  end function needCalcQValByBCID
  ! ****************************************************************************!

  ! ************************************************************************** !
  !> Obtain the surface normal of the chosen boundary if the boundary has
  !! `store_normal` set. Otherwise a null vector is returned.
  function getNormal(geometry, elemBary, node, node_pos, bcid) result(normal)
    ! --------------------------------------------------------------------- !
    type(sdr_geometry_type), intent(in) :: geometry
    real(kind=rk), intent(in) :: elemBary(3)
    type(sdr_node_type), intent(in) :: node
    integer, intent(in) :: node_pos
    integer, intent(in) :: bcid
    real(kind=rk) :: normal(3)
    ! --------------------------------------------------------------------- !
    logical :: calc_normal
    integer :: iAtt
    integer :: intersected_first, intersected_last
    integer :: geom_prim, atb_pos, prim_pos
    integer :: obj_pos
    integer :: iSpatialObj
    integer :: iTriangle
    integer :: nIntersected
    integer :: nAdjacent
    real(kind=rk) :: tridistance
    type(grw_realArray_type) :: closedist
    integer, allocatable :: closekind(:)
    integer, allocatable :: distrank(:)
    real(kind=rk), allocatable :: closepoint(:,:)
    real(kind=rk), allocatable :: trinormal(:,:)
    real(kind=rk) :: mindist
    ! --------------------------------------------------------------------- !

    normal = 0.0_rk
    calc_normal = .false.

    if ((bcid <= geometry%attribute%kindpos(sdr_Boundary_object)%nVals) &
      & .and. (bcid >= 1)) then
      iAtt = geometry%attribute%kindpos(sdr_Boundary_object)%val(bcid)
      calc_normal = geometry%attribute%dynArray%val(iAtt)%store_normal
    end if

    if (calc_normal) then
      intersected_first = node%userObjPos%val(node_pos)%first
      intersected_last = node%userObjPos%val(node_pos)%last
      nIntersected = intersected_last - intersected_first + 1
      call init(closedist, length = nIntersected)
      allocate(closekind(nIntersected))
      allocate(closepoint(3, nIntersected))
      allocate(trinormal(3, nIntersected))
      do iSpatialObj = intersected_first, intersected_last
        obj_pos = node%intersected_object%val(iSpatialObj)
        atb_pos = geometry%spatialObj%val(obj_pos)%attribute_position
        same_attr: if (atb_pos == iAtt) then
          geom_prim = geometry%spatialObj%val(obj_pos)%geometry_primitive
          prim_pos = geometry%spatialObj%val(obj_pos)%primitive_position
          ! TODO: find closest point on the geometrical objects to the
          !       element bary center.
          select case(geom_prim)
          case(triangle)
            iTriangle = closedist%nVals + 1
            call tem_triangle_normal_proximity(          &
              &    me = geometry%triangle%val(prim_pos), &
              &    point = elemBary,                     &
              &    closest = closePoint(:,iTriangle),    &
              &    closekind = closekind(iTriangle),     &
              &    distance = tridistance,               &
              &    normal   = trinormal(:,iTriangle)     )
            call append(closedist, val = tridistance)
          end select
        end if same_attr
      end do

      if (closedist%nVals > 0) then
        normal = trinormal(:,1)
        if (closedist%nVals > 1) then
          allocate(distrank(closedist%nVals))
          ! Rank the closest points of all triangles by their distance
          call mrgrnk(closedist%val(:closedist%nVals), distrank)
          normal = trinormal(:,distrank(1))

          if (closekind(distrank(1)) > tem_triangle_close_face) then
            ! The closest point on the triangle lies on the boundary of the
            ! triangle. We also need to consider adjacent triangles that may
            ! share this point and average their normals.
            mindist = closedist%val(distrank(1)) * (1.0_rk + epsilon(mindist))
            nAdjacent = 1
            do iTriangle=2,closedist%nVals
              if (closedist%val(distrank(iTriangle)) > mindist) EXIT
              normal = normal + trinormal(:,distrank(iTriangle))
              nAdjacent = nAdjacent + 1
            end do
            normal = normal / nAdjacent
            normal = dot_product(elemBary - closePoint(:,distrank(1)), normal) * normal
          end if
          deallocate(distrank)

        end if
      end if

      deallocate(closekind)
      deallocate(closepoint)
      deallocate(trinormal)
      call destroy(closedist)

    end if

  end function getNormal
  ! ************************************************************************** !

  ! ****************************************************************************!
  !> This routine checks if a boundary need flood periphery for diagonal
  !! directions for a given BCID.
  !! It is used in identify_boundary routine
  function needFldDglByBCID( attribute, bcid ) result( flood_diagonal )
    type(sdr_attrList_type), intent(in) :: attribute
    integer, intent(in) :: bcid
    integer :: iAtt
    logical :: flood_diagonal

    if ((bcid <= attribute%kindpos(sdr_Boundary_object)%nVals) &
      & .and. (bcid >= 1)) then
      iAtt = attribute%kindpos(sdr_Boundary_object)%val(bcid)
      flood_diagonal = attribute%dynArray%val(iAtt)%flood_diagonal
    else
      flood_diagonal = .false.
    end if
  end function needFldDglByBCID
  ! ****************************************************************************!

  ! ****************************************************************************!
  !> This routine gives special treatment when qVal > 1.0 or qVal == -1.0
  !! for flooded neighbor, treat it as normal fluid: clean BCID,
  !! set qVal to -1 (no itersection).
  !! for non-flooded neighbor, treat it as high order wall: set qVal to 1
  subroutine sdr_truncate_qVal( proto, qVal, BCID, neighPos, treeID_periodic )
    ! --------------------------------------------------------------------------!
    !> preliminary tree
    type(sdr_protoTree_type), intent(in) :: proto
    !> qValue
    real(kind=rk),intent(inout) :: qVal
    !> boundary id
    integer(kind=long_k), intent(inout) :: BCID
    !> neighbor position in proto tree
    integer, intent(in) :: neighPos
    !> negative treeID of periodic domain
    integer(kind=long_k), intent(in), optional :: treeID_periodic
    ! --------------------------------------------------------------------------!

    if( qVal > 1._rk .or. qVal == sdr_qVal_no_intersection ) then
      if (sdr_nodeProp_btest(node  = proto%node,  &
        &                    iNode = neighPos,    &
        &                    bit   = isFlooded_bit) ) then
        ! clear BCID or set periodic boundary treeID
        if (present(treeID_periodic)) then
          BCID = - treeID_periodic
        else
          BCID = 0_long_k
        endif
        ! set qVal to no intersection
        qVal = sdr_qVal_no_intersection
      else
        qVal = 1._rk
      end if
    end if
  end subroutine sdr_truncate_qVal
  ! ****************************************************************************!

end module sdr_boundary_module
