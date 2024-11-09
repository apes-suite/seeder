! Copyright (c) 2012-2014, 2016 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2012-2015 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2012, 2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2012, 2014 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
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
! *****************************************************************************
!> This module provides the functionality to find the part of the universe
!! cube, which is actually building up the computational domain.
module sdr_flooding_module
  use env_module,           only: rk, long_k, eps
  use tem_param_module,     only: qQQQ, qOffset
  use treelmesh_module,     only: treelmesh_type
  use tem_geometry_module,  only: tem_eligibleChildren, tem_ElemSize,          &
    &                             tem_baryOfID
  use tem_topology_module,  only: tem_directChildren, tem_LevelOf,             &
    &                             tem_parentOf, tem_CoordOfID, tem_IdOfCoord
  use tem_tools_module,     only: tem_positionInSorted
  use tem_logging_module,   only: tem_log, tem_toStr
  use tem_debug_module,     only: main_debug

  use sdr_protoTree_module, only: sdr_protoTree_type, sdr_neighbor_in_proto,   &
    &                             sdr_write_proto_as_restart
  use sdr_geometry_module,  only: sdr_geometry_type
  use sdr_node_module,      only: isFlooded_bit, IntersectsBoundary_bit,    &
    &                             sdr_wetNeighborsFace, sdr_mark_floodNode, &
    &                             isLeaf_bit
  use sdr_config_module,    only: sdr_confHead_type
  use sdr_boundary_module,  only: needCalcQValByBCID, sdr_qValByNode,          &
    &                             needFldDglByBCID, sdr_qVal_no_intersection

  use sdr_timer_module,     only: timer_handle_flooding
  use tem_timer_module,     only: tem_startTimer,tem_stopTimer 


  implicit none

  private

  public :: sdr_flood_tree

contains

  ! *****************************************************************************
  !> This routine identifies the nodes, which are supposed to be part of the
  !! computational domain, as defined by the seed objects.
  !!
  subroutine sdr_flood_tree( proto, geometry, header, meshUniverse )
    ! --------------------------------------------------------------------------!
    !> The proto tree description with all the data enabling the flooding.
    type(sdr_protoTree_type), intent(inout) :: proto
    !> Description of geometric objects. Propably not needed here, remove again
    !! if this is the case.
    type(sdr_geometry_type), intent(in) :: geometry
    !> some global information on solver name and version
    type(sdr_confHead_type), intent(inout) :: header
    !> treelmesh contains bounding cube info
    type(treelmesh_type), intent(in) :: meshUniverse
    ! --------------------------------------------------------------------------!
    ! --------------------------------------------------------------------------!

    call tem_startTimer( timerHandle = timer_handle_flooding ) 

    call tem_log(1, 'Flooding tree ... ')

    call floodwaves_tree( proto    = proto,    &
      &                   header   = header,   &
      &                   geometry = geometry  )

    call tem_log(2, 'Number of Flooded leaves: ' &
      &             // trim(tem_toStr(proto%nFloodedLeaves)) )

    ! JQ: if qVal is required, flood some additional nodes. only do it once.
    if ( any(geometry%attribute%dynArray%val( : )%calc_dist) ) then
      call flood_periphery( proto, geometry, meshUniverse )
    end if

    !! Mark all virtual nodes, which contain a flooded child as flooded
    !! starting from the second finest level moving up to the root.
    !! This allows to easily avoid non-flooded domains later on.
    call flood_parents( proto, geometry%attribute%color_inverted )

    call tem_stopTimer( timerHandle = timer_handle_flooding )

  end subroutine sdr_flood_tree
  ! *****************************************************************************


  ! *****************************************************************************
  !> This routine loop over all nodes are flood non-interesting leaf node with
  !! wet face and inherit the wetness of the virtual node to the eligble
  !! childrens
  !!
  !! The algorithm works by flooding the domain, starting from the seed points
  !! up to boundary elements. This approach is quite robust to broken STL
  !! definitions, as any cracks below the resolution are automatically healed
  !! and there is no dependece on the orientation of the surfaces.
  !! To avoid unintended leaking, the flooding takes only the 6 side neighbors
  !! into account, that is the computational domain will always be connected
  !! by faces, there will be no parts of the domain which are only connected
  !! by edges or corners.
  !!
  !! Seeds are already marked as flooded during the previous leaf
  !! identification, also the neighboring sides have already been marked as wet.
  !! Note, that seeds in nodes with boundaries are ignored and not flooded.
  !!
  !! Several iterations are done, referred to as waves and within each the
  !! following algorithm is used:
  !!
  !! Iterate over all nodes.
  !! - If node is a leaf
  !!   + If not intersecting boundaries:
  !!     * Check if any side is wet, and if so, flood yourself.
  !!       If flooded, mark all neighboring sides as wet.
  !!   + If intersecting boundaries: nothing to do!
  !!
  !! - For virtual nodes: (if not leaf)
  !!   + Inherit wet sides down to the direct children on this side
  !!     (eligible children).
  !!
  !! The procedure is finished, if no property changed during a wave.
  !! (Flooded status and wet faces are encoded in the PropertyBits field.
  !!
  subroutine floodwaves_tree( proto, header, geometry )
    ! --------------------------------------------------------------------------!
    !> The proto tree description with all the data enabling the flooding.
    type(sdr_protoTree_type), intent(inout) :: proto
    !> some global information on solver name and version
    type(sdr_confHead_type), optional, intent(in) :: header
    type(sdr_geometry_type), intent(in) :: geometry
    ! --------------------------------------------------------------------------!
    integer :: iNode, faceBitCheck, iDir, iSide
    integer(kind=long_k) :: childIds(8)
    integer(kind=long_k) :: virtualId
    integer, allocatable :: eligible_childs(:)
    integer :: childPos, wetFace
    integer :: wetFaceBit_pos
    integer :: iChild
    integer, allocatable :: old_prop(:,:)
    integer :: iWave
    integer :: iColor
    integer :: col_int
    integer :: col_bit
    integer :: myColor
    integer :: nodeprops
    integer :: propLength
    ! --------------------------------------------------------------------------!

    propLength = proto%node%propLength
    allocate(old_prop(propLength, proto%node%nNodes))

    ! The layout of each color is:
    ! intersection with boundaries is set in bit 0
    ! volume flooded is indicated by bit 1
    ! faces are indicated by bits 2-7
    ! Thus, if any face is wet, the color is >= 4
    faceBitCheck = 4

    ! Iterations to flood (waves)
    ! (limit by number of nodes, we should not need to do more than number
    !  of nodes iterations, though this limit is somewhat arbitrary)
    ! Leave loop, when no property was changed during the wave.
    ! This should typically happen much faster.
    waves: do iWave=1,proto%node%nNodes
      ! First save the current properties for later comparison:
      old_prop = proto%node%PropertyBits%val

      ! Iterate over all nodes.
      ! Nodes are sorted level-wise with the coarser levels coming first.
      ! The nodeloop thus iterates over the elements in a top down fashion,
      ! wet faceness of coarser nodes is therefore automatically propagated down
      ! within a single wave.
      nodeLoop: do iNode = 1, proto%node%nNodes

        nodeprops = ibits(proto%node%PropertyBits%val(propLength, iNode), &
          &               pos = proto%node%lastbyte_pos, len=8)

        ! Depending on the type of node take some action
        isleaf: if ( btest(nodeprops, isLeaf_Bit) ) then
          ! For leaves:

          col_loop1: do iColor=1,proto%node%nColors
            col_int = (iColor-1) / proto%node%bytes_per_int + 1
            col_bit = mod(iColor-1, proto%node%bytes_per_int)*8
            myColor = ibits(proto%node%PropertyBits%val(col_int, iNode), &
              &             pos = col_bit, len = 8)

            ! myColor bit 0 encodes the intersecting status of the node
            ! myColor bit 1 encodes its flooding, both have to be not true,
            ! otherwise we have nothing to do. Thus, the following check
            ! ensures, they are both 0:
            if (mod(myColor,4) == 0) then
              ! NOT intersecting (relevant) boundaries, and not flooded yet.

              ! Check if any of my sides are wet, and if so, flood myself and
              ! wet my neighbors sides accordingly.
              if (myColor >= faceBitCheck) then

                ! I got a wet side and am a leaf node, flood myself now and
                ! increase the flooded leaf counter accordingly.
                call sdr_mark_floodNode( proto%node,           &
                  &                      iNode,                &
                  &                      proto%nFloodedLeaves, &
                  &                      color = iColor )

                ! Now mark all neighboring sides with the color set above.
                call sdr_wetNeighborsFace(proto%node, iNode, col_int, col_bit)

              end if ! at least one of my faces is wet

            end if

          end do col_loop1

        else isleaf

          ! For virtual nodes:
          ! Inherit wet sides down to the direct children on this side
          ! (eligible children).
          ! Only need to take action, if any of the faces is wet.
          col_loop2: do iColor=1,proto%node%nColors
            col_int = (iColor-1) / proto%node%bytes_per_int + 1
            col_bit = mod(iColor-1, proto%node%bytes_per_int)*8
            myColor = ibits(proto%node%PropertyBits%val(col_int,iNode), &
              &             pos = col_bit, len = 8)
            wetside: if (myColor >= faceBitCheck) then
              ! Find direct childern of the virtual node and then check for
              ! specific wet face and find the eligible children in the
              ! direction of the wet face and mark the corresponding faces in
              ! the neighbor in that direction as wet.

              ! Some auxilary variables
              virtualId = proto%node%TreeID%val( iNode )
              childIds = tem_directChildren( virtualId )

              ! Iterate over all 6 sides of the node
              do iDir = 1,3
                do iSide = 1,2

                  wetFace = ( (iDir-1)*2 + iSide ) + 1
                  wetFaceBit_pos = col_bit + wetFace

                  ! Check if this face is wet
                  if (btest(myColor, wetFace))  then
                    ! Encountered a wet side, select the 4 elegible children
                    ! adjacent to this side.
                    ! (Have to use the treelm numbering of sides, given by the
                    !  iDir + 3*(iSide-1) formula).
                    call tem_eligibleChildren( eligible_childs, &
                      &                        iDir + 3*(iSide-1) )

                    ! Iterate over all children and mark their faces in the
                    ! current direction as wet.
                    do iChild = 1, size(eligible_childs)
                      ! The first child is given by linkpos(1), and all others
                      ! follow consecutively.
                      childpos = proto%node%linkpos(1)%val(iNode) &
                        &      + eligible_childs(iChild) - 1

                      ! Wet this childs face in the current direction.
                      proto%node%PropertyBits%val(col_int, childPos) &
                        & = ibset(proto%node%PropertyBits &
                        &                   %val(col_int, childPos), &
                        &         wetFaceBit_pos)
                    end do
                  end if

                end do ! side loop
              end do ! dirloop

            end if wetside

          end do col_loop2

        end if isleaf

      end do nodeLoop

      ! output debug information as restart files
      if( (trim(main_debug%debugMesh) .ne. '') .and. &
        & (main_debug%debugFiles) ) then
        !debug output
        call sdr_write_proto_as_restart(proto, geometry, &
          &                             iwave, header, 'flood')
      end if

      ! Compare the new properties to the old ones, if nothing changed leave the
      ! loop!
      if (all(ieor(old_prop, proto%node%PropertyBits%val) == 0)) EXIT waves

    end do waves

    deallocate(old_prop)

    call tem_log(1, '               ... done with nWaves = ' &
      &             // trim(tem_toStr(iWave)) )

  end subroutine floodwaves_tree
  ! *****************************************************************************


  ! *****************************************************************************
  !> This routine loops over all intersected with geoemtry nodes and fluidify
  !! some node according to the following rule:\n
  !! 1. one of its link does noe intersect with any geometry that requires qVal
  !! 2. it has fluid neighbor on that direction.
  !!    i.e. it is wet in that side.
  !! Jiaxing Qi
  !!
  !!@todo HK: works for single color only right now!
  !!          (qvalues are only computed for first color, need to think about
  !!           what to do for multiple colors.)
  subroutine flood_periphery( proto, geometry, meshUniverse )
    ! --------------------------------------------------------------------------!
    !> The proto tree description with all the data enabling the flooding.
    type(sdr_protoTree_type), intent(inout) :: proto
    !> Description of geometric objects.
    type(sdr_geometry_type), intent(in) :: geometry
    !> treelmesh contains bounding cube info
    type(treelmesh_type), intent(in) :: meshUniverse
    ! --------------------------------------------------------------------------!
    integer :: iNode, iDir, iSide, iLink, wetFace
    integer :: wetFaceBit_pos
    integer(kind=long_k) :: treeID
    !> qVal for 6 faces
    real(kind=rk) :: qVal( 1:qQQQ )
    real(kind=rk) :: bary(3)
    real(kind=rk) :: dx
    integer :: oldFlood
    integer :: myColor
    integer :: col_int, col_bit
    integer :: iColor
    integer :: nodeprops
    integer :: proplength
    ! --------------------------------------------------------------------------!

    proplength = proto%node%propLength
    oldFlood = proto%nFloodedLeaves
    call tem_log(1, 'Flooding periphery (calc_dist active) ... ')

    ! loop over all nodes
    do iNode = 1, proto%node%nNodes
      nodeprops = ibits(proto%node%PropertyBits%val(propLength,iNode), &
        &               pos = proto%node%lastbyte_pos, len=8)
      ! Only proceed if this is a leaf.
      isLeaf: if ( btest(nodeprops, isLeaf_Bit) ) then

        ! Act on qValues only in the first color!
        do iColor=1,1 !proto%node%nColors
          col_int = (iColor-1) / proto%node%bytes_per_int + 1
          col_bit = mod(iColor-1, proto%node%bytes_per_int)*8
          myColor = ibits(proto%node%PropertyBits%val(col_int,iNode), &
            &             pos = col_bit, len = 8)

          ! If the node is not flooded yet (bit 1) and intersects a boundary
          ! relevant to this color (bit 0)
          nonfloodbc: if (mod(myColor,4) == 1) then

            ! Now check the intersected boundaries for the need to find the
            ! q-values.
            needq: if ( needCalcQValByBCID(geometry%attribute, &
              &                     proto%node%minBCID%val(iNode)) ) then

              treeID = proto%node%treeID%val( iNode )
              Bary = tem_BaryOfID( meshUniverse, treeID )
              dx = tem_ElemSize( meshUniverse, treeID )
              qVal = sdr_qVal_no_intersection
              ! loop over all 26 directions as per order in tem_param_module
              ! compute all qVal once and use it for axis normal and
              ! diagonal axis
              do iDir = 1, qQQQ
                  ! calculate qVal from myself (boundarz node)
                  call sdr_qValByNode( proto, geometry, dx, iDir, &
                    & Bary, iNode, qVal(iDir) )
              end do

              ! Only proceed if no boundary runs through the bary center
              ! (qval == 0)
              qValEpsCheck: if (minval(qVal, mask=qVal>=0) > eps) then

                ilinkLoop: do iDir = 1,3
                  do iSide = 1,2

                    ! loop over 6 flooding links: -x, +x, -y, +y, -z, +z
                    iLink = iDir + 3*(iSide-1)
                    wetFace = iSide + (iDir-1)*2
                    wetFaceBit_pos = wetface + col_bit + 1

                    ! When there is no intersection in this direction, we need
                    ! to check if flooding of the current node is necessary.
                    if (qVal(iLink) < 0.0_rk) then
                      ! If the neighbor in this direction is a fluid, it has
                      ! wetted my side. As it is reachable without any obstacle
                      ! in between, the current node needs to be flooded, too.
                      if ( btest(proto%node%PropertyBits%val(col_int,iNode), &
                        &        wetfaceBit_pos) ) then
                        call sdr_mark_floodNode( proto%node,           &
                          &                      iNode,                &
                          &                      proto%nFloodedLeaves, &
                          &                      color = iColor        )

                        exit iLinkLoop
                      end if
                    end if
                  end do
                end do ilinkLoop

                ! Check qval for diagonal directions for flooding
                if ( needFldDglByBCID(geometry%attribute, &
                  &                   proto%node%minBCID%val(iNode)) ) then
                  call flood_periphery_diagonal( proto, iNode, treeID, &
                    &                            qVal, iColor          )
                end if

              end if qValEpsCheck

            end if needq

          end if nonfloodbc

        end do

      end if isLeaf
    end do ! iNode

    call tem_log(1, ' ... Done flooding the periphery!')

    call tem_log(2, 'Number of Flooded leaves in periphery: '           &
      &             // trim(tem_toStr(proto%nFloodedLeaves - oldFlood)) )

    call tem_log(1, 'Final number of Flooded leaves: '       &
      &             // trim(tem_toStr(proto%nFloodedLeaves)) )

  end subroutine flood_periphery
  ! *****************************************************************************


  ! *****************************************************************************
  !> This routine checks for qVal of the periphery and floods if qVal < 0 and
  !! the node in that direction is fluid and not intersected by boundary
  subroutine flood_periphery_diagonal( proto, node_pos, treeID, qVal, iColor )
    ! --------------------------------------------------------------------------!
    !> The proto tree description with all the data enabling the flooding.
    type(sdr_protoTree_type), intent(inout) :: proto
    !> node position in protoTree
    integer, intent(in) :: node_pos
    !> treeID of current node
    integer(kind=long_k), intent(in) :: treeID
    !> qVal for all 26 neighbor directions, should be calculated already
    real(kind=rk) :: qVal( : )
    !> Color to do the flooding in
    integer, intent(in) :: iColor
    ! --------------------------------------------------------------------------!
    integer :: coord(4), neighbor_pos, neighbor_level, iDir
    integer :: col_int, col_bit, myColor
    integer :: neighColor
    ! --------------------------------------------------------------------------!

    col_int = (iColor-1) / proto%node%bytes_per_int + 1
    col_bit = mod(iColor-1, proto%node%bytes_per_int)*8
    myColor = ibits(proto%node%PropertyBits%val(col_int, node_pos), &
      &             pos = col_bit, len = 8)

    coord = tem_coordOfId( treeID = treeID )
    off_dirLoop: do iDir = 7, qQQQ

      ! fluidify this node, if not intersected with boundary i.e. qVal = -1.0
      ! and has a fluid in that direction
      if (qVal(iDir) < 0) then

        ! get neighbor position in protoTree
        neighbor_pos = sdr_neighbor_in_proto( proto, coord, iDir, &
          &                                   neighbor_level )

        neighColor = ibits(proto%node%PropertyBits%val(col_int, neighbor_pos), &
          &                pos = col_bit, len = 8)

        ! check if neighbor in this direction is flooded (bit 1) and
        ! not intersected by any geometry (bit 0)
        if ( mod(neighColor,4) == 2) then

          ! mark it as isFlooded
          ! these nodes are flooded then, but also intersect boundaries
          ! Found a new flooded leaf, increase the counter accordingly.
          call sdr_mark_floodNode( node           = proto%node,           &
            &                      iNode          = node_pos,             &
            &                      nFloodedLeaves = proto%nFloodedLeaves, &
            &                      color          = iColor                )

          return
        end if
      end if
    end do off_dirLoop

  end subroutine flood_periphery_diagonal
  ! *****************************************************************************


  ! *****************************************************************************
  !> Mark all virtual nodes, which contain a flooded child as flooded
  !! starting from the second finest level moving up to the root.
  !! This allows to easily avoid non-flooded domains later on.
  subroutine flood_parents( proto, color_inverted )
    ! --------------------------------------------------------------------------!
    !> The proto tree description with all the data enabling the flooding.
    type(sdr_protoTree_type), intent(inout) :: proto
    !> List of flags for each color to indicate, wether the color should be
    !! inverted after flooding.
    logical, intent(in) :: color_inverted(:)
    ! --------------------------------------------------------------------------!
    integer :: iNode, iLevel, firstchild
    integer :: childprops(8)
    integer :: nodeflooded_bit
    integer :: propLength
    integer :: nodeprops
    integer :: col_int, col_bit
    integer :: iColor
    integer :: nInverted
    integer :: iInv
    integer, allocatable :: invColor(:)
    ! --------------------------------------------------------------------------!
    propLength = proto%node%propLength
    nodeflooded_bit = proto%node%lastbyte_pos + isFlooded_bit

    nInverted = count(color_inverted)
    allocate(invColor(nInverted))
    ! If any color needs to be inverted, run through the finest level also, and
    ! invert color floodings as needed.
    if (nInverted>0) then
      iInv = 0
      do iColor=1,proto%node%nColors
        if (color_inverted(iColor)) then
          iInv = iInv + 1
          invColor(iInv) = iColor
        end if
      end do
      iLevel = proto%nLevels
      ! All nodes on the finest level have to be leaves.
      do iNode = proto%levelNode_first(iLevel), proto%levelNode_last(iLevel)
        nodeprops = ibits(proto%node%PropertyBits%val(propLength,iNode), &
          &               pos = proto%node%lastbyte_pos, len=8)

        ! Only consider elements, that are flooded by any color (especially the
        ! none color). If an element is not belonging to the domain identified
        ! positively by the none color, we will not invert its color floodings.
        if ( btest(nodeprops, isFlooded_bit) ) then
          do iInv=1,nInverted
            col_int = (invColor(iInv)-1) / proto%node%bytes_per_int + 1
            col_bit = mod(invColor(iInv)-1, proto%node%bytes_per_int)*8
            if ( btest(proto%node%PropertyBits%val(col_int, iNode), &
                &     col_bit+1) ) then
              proto%node%PropertyBits%val(col_int, iNode)               &
                &  = ibclr(proto%node%PropertyBits%val(col_int, iNode), &
                &          col_bit+1)
            else
              proto%node%PropertyBits%val(col_int, iNode)               &
                &  = ibset(proto%node%PropertyBits%val(col_int, iNode), &
                &          col_bit+1)
            end if
          end do
        end if
      end do
    end if

    do iLevel = proto%nLevels-1, 0, -1
      do iNode = proto%levelNode_first(iLevel), proto%levelNode_last(iLevel)
        nodeprops = ibits(proto%node%PropertyBits%val(propLength,iNode), &
          &               pos = proto%node%lastbyte_pos, len=8)

        if ( btest(nodeprops, isLeaf_bit) ) then
          ! Check leaves for inversion of flooding status.
          if ( btest(nodeprops, isFlooded_bit) ) then
            do iInv=1,nInverted
              col_int = (invColor(iInv)-1) / proto%node%bytes_per_int + 1
              col_bit = mod(invColor(iInv)-1, proto%node%bytes_per_int)*8
              if ( btest(proto%node%PropertyBits%val(col_int, iNode), &
                  &     col_bit+1) ) then
                proto%node%PropertyBits%val(col_int, iNode)               &
                  &  = ibclr(proto%node%PropertyBits%val(col_int, iNode), &
                  &          col_bit+1)
              else
                proto%node%PropertyBits%val(col_int, iNode)               &
                  &  = ibset(proto%node%PropertyBits%val(col_int, iNode), &
                  &          col_bit+1)
              end if
            end do
          end if

        else
          ! If any of the node children are flooded, the node has to be
          ! considered in the final mesh. Thus, mark it accordingly in this
          ! case. Only non-leaves need to be considered here.
          firstChild = proto%node%linkPos(1)%val(iNode)

          childprops = ibits(proto%node%PropertyBits &
            &                     %val(propLength, firstChild &
            &                                      :firstChild+7), &
            &                pos = proto%node%lastbyte_pos, len=8)
          if ( any(btest(childprops, isFlooded_bit)) ) then
            proto%node%PropertyBits%val(propLength,iNode) &
              &  = ibset(proto%node%PropertyBits%val(propLength, iNode), &
              &          nodeflooded_bit)
          end if

          do iColor=1,proto%node%nColors
            col_int = (iColor-1) / proto%node%bytes_per_int + 1
            col_bit = mod(iColor-1, proto%node%bytes_per_int)*8
            childprops = ibits( proto%node%PropertyBits%val(col_int,        &
              &                                             firstChild      &
              &                                             :firstChild+7), &
              &                 pos = col_bit, len = 8                      )
            if ( any(btest(childprops, 1)) ) then
              proto%node%PropertyBits%val(col_int,iNode)                &
                &  = ibset(proto%node%PropertyBits%val(col_int, iNode), &
                &          col_bit+1)
            end if
          end do

        end if

      end do
    end do

  end subroutine flood_parents
  ! *****************************************************************************

 
end module sdr_flooding_module
