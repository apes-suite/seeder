! Copyright (c) 2011-2016 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2011 Daniel Harlacher <d.harlacher@grs-sim.de>
! Copyright (c) 2011, 2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2011-2012 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2011-2015 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2011-2012 Metin Cakircali <m.cakircali@grs-sim.de>
! Copyright (c) 2011 Aravindh Krishnamoorthy <aravindh28.4@gmail.com>
! Copyright (c) 2011 Khaled Ibrahim <k.ibrahim@grs-sim.de>
! Copyright (c) 2011-2012 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2012 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2015 Samuel Ziegenberg
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2016 Raphael Haupt <Raphael.Haupt@student.uni-siegen.de>
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
!> S E E D E R
!! Mesh generator within the APES suite
!!
!! For a documentation, run ./waf gendoxy and find the documentation at
!! ./Documentation/html/index.html
program seeder
  ! treelm modules
  use tem_general_module,      only: tem_start, tem_finalize
  use tem_debug_module,        only: main_debug
  use tem_timer_module,        only: tem_stopTimer

  ! seeder modules
  use sdr_aux_module,          only: sdr_init_global
  use sdr_config_module,       only: sdr_load_config, sdr_confHead_type
  use sdr_geometry_module,     only: sdr_geometry_type
  use sdr_protoTree_module,    only: sdr_protoTree_type, sdr_build_protoTree,  &
    &                                sdr_write_proto_as_restart
  use sdr_flooding_module,     only: sdr_flood_tree
  use sdr_proto2treelm_module, only: sdr_temData_type, sdr_proto2treelm,       &
    &                                sdr_dump_treelm
  use sdr_refinePT_module,     only: sdr_refine_leaf, sdr_smooth_leaf,         &
    &                                sdr_inHerit_distanceRefineObject

  use sdr_timer_module,        only: sdr_addTimers, sdr_dumptimers

  implicit none

  ! -----------------------------------------------------------------------------
  type(sdr_confHead_type) :: header !< global header variable
  type(sdr_geometry_type) :: geometry !< contains all geometry object
  type(sdr_protoTree_type) :: protoTree !< preliminary tree
  type(sdr_temData_type) :: temData !< TreElM data to produce
  ! -----------------------------------------------------------------------------


  ! Initialize environment.
  call tem_start(codeName = 'Seeder',      &
    &            general  = header%general )
  
  ! add a timer object to measure time in different steps of seeder 
  call sdr_addTimers()

  if (header%general%proc%rank == 0) then
   ! Initialize seeder screen and solveHead type.
   call sdr_init_global(header%general%solver)
  end if

  ! Load configuration data according to command line arguments.
  call sdr_load_config(header, geometry)

  ! Build the preliminary prototree to narrow down the boundaries of the
  ! computational domain down to the desired resolution.
  call sdr_build_protoTree(protoTree, geometry, header)

  !copy bounding cube info to temData
  temData%meshUniverse%global%origin = geometry%universe%origin
  temData%meshUniverse%global%BoundingCubeLength = geometry%universe%extent

  if (main_debug%active) &
    &  call sdr_write_proto_as_restart(protoTree, geometry, 0, header )

  ! Flood the prototree, starting from the seed, to identify all elements that
  ! are to be part of the computational domain.
  call sdr_flood_tree( protoTree, geometry, header, temData%meshUniverse )

  ! Dump proto tree if requested.
  if (main_debug%active) &
    &  call sdr_write_proto_as_restart(protoTree, geometry, 1, header, &
    &                                  'flooded')

  ! Inherit intersected distance refine objects
  if (protoTree%node%distanceRefine) &
    & call sdr_inHerit_distanceRefineObject( protoTree, geometry )

  ! Refine protoTree leafs to its finest level defined by refinement
  ! object or minlevel
  call sdr_refine_leaf( protoTree, geometry )

  ! Dump proto tree if requested.

  ! Smooth protoTree to avoid level jump > 1
  if (geometry%smoothLevels) then
    ! Dump protoTree before smooth setup if debug is active
    if (main_debug%active) &
      &  call sdr_write_proto_as_restart(protoTree, geometry, 2, header, &
      &                                  'refined')

    call sdr_smooth_leaf(protoTree, header, protoTree%nLevels) 
  end if  

  ! Dump proto tree if requested.
  if (main_debug%active) &
    &  call sdr_write_proto_as_restart(protoTree, geometry, 100, header, &
    &                                  'final')

  ! Everything is ready now to define the computational mesh itself.
  ! For this we have to refine within defined refinement objects and
  ! identify boundary conditions.
  call sdr_proto2treelm(protoTree, geometry, temData, header)

  ! Write the treelmesh and boundary list to disk.
  call sdr_dump_treelm(temData, geometry, header)

  ! stop seeder timer before dump timing
  call tem_stopTimer(timerHandle = header%general%solver%timerHandle)

  ! dump timers to timing output file
  call sdr_dumptimers(general = header%general,        &
    &                 nFluids = temData%treeID%nVals,  &
    &                 nBnds   = temData%BC_ID(1)%nVals )

  ! Finalize environment.
  call tem_finalize(header%general)

end program seeder


!>\page sdr_algorithm Seeder Algorithm
!!\author Kannan Masilamani
!!
!! In a first step, Seeder will create a protoTree, that refines all domain
!! boundaries with the required resolution. This is done iteratively level by
!! level, and each node, that is intersected by a domain boundary is further
!! refined on the next level, if the maximum level of all intersected objects
!! or global minimum level is not reached yet.
!! At the end of this process the boundaries of the domain are resolved as
!! requested everywhere.
!!
!! After the boundaries are properly resolved, it needs to be decided, what is
!! actually part of the domain itself. For this elements are marked as flooded
!! based upon their neighborhood. Obviously this only works properly if there
!! is some initial starting point to get the flooding running. This is achieved
!! by the user defined seed objects, which will mark all the elements that they
!! intersect as flooded initially, as long as the elements do not intersect
!! boundaries at the same time. Elements which intersect boundaries are never
!! flooded!.
!!
!! Before flooding, the neighbors of each leaf node is identified for 
!! each node faces which is necassary for flooding to find the state of 
!! neighbor node.
!!
!! This flooding is done iteratively until a iteration did not change any
!! element state anymore. All flooded elements identify the computational
!! domain now.
!!
!! Finally the computational domain is refined to the desired level everywhere,
!! as given by refinement objects or the global minimal refinement level and
!! boundary conditions are assigned.
!! This data is then written in treelm format to disk at the location, specified
!! by the user in the configuration file. Note, that if this is a directory,
!! you should not forget the trailing path separator in the definition and
!! create the directory beforehand.

!>\page Boundary conditions
!!
!! Boundary conditions are defined by the label attached to a boundary object.
!! Multiple boundary objects might have the same label and therefore attached
!! to the same boundary condition.
!! If there are multiple boundaries in a given direction of an element, the
!! one defined first in the configuration, will be used. This allows you to
!! set the precedence order of boundaries that should be used at any
!! intersections.
!!
!! Only boundary labels that actually appear in the final mesh will be written
!! to the header description of the boundaries in the mesh.

!>\page Colored seeds
!!
!! Seeds might have colors and thereby define different areas within the mesh.
!! This is for example useful to ascribe special material parameters to parts
!! of the domain.
!! Each colored seed is bounded by a boundary of the same color, or boundaries
!! which have the color 'none' or 'all'.
!! Elements which have a non 'none' color attached to them will have the
!! isColored property and the additional property information about which colors
!! are attached to that element will be stored.
!! There might be an arbitrary number of colors and multiple colors can be
!! assigned to each element.
!! Color names can be chosen arbitrarily and are case insensitive.
!!
!! We achieve this by using a bitfield for each color and basically do the
!! flooding for each color. With up to three colors there should be little to
!! no memory overhead imposed by this approach, as each color requires only one
!! byte in the integer bitfield, that we use (the last byte is reserved for
!! general node information). Beyond that additional integers will be used as
!! needed by the number of colors.
!! In the mesh on the disk, colors are stored in ASCII characters with 7 colors
!! per byte, and meshes without coloring ('none' color everywhere as only color)
!! no color information will be written at all.
!!
!! Note: colored boundaries are subject to the same rules of boundary labels
!! as none colored boundaries and if they intersect with none colored boundaries
!! will take precedence over those, if configured first.
!! As this is typically not desired, you probably want to define your colored
!! boundary objects after the none colored ones!

!>\page Subresolved boundaries
!!
!! Boundaries with a color might be resolved beyond the stair-case representation
!! offered by the mesh by creating a polynomial representation of the color
!! values within the intersected elements.
!! For this you need to define a subresolution table, describing at least the
!! polynomial degree that should be used for this information.
!! When this table is defined, you can indicate that a (colored) boundary is to
!! be subresolved by stating 'subresolution = true' in its attribute.
!!
!! Polynomial information will be created per color and is accessible as a
!! property in the mesh afterwards.

!>\page Distance refinement
!!
!! 

!>\page Smoothing

!> \page special Special geometries
!! Seeder can create some special geometries. 
!! Currently is supports only one special geometry, spacer.
!! \li \ref spacers "Spacers"
!!
!> \example testsuite/plane/seeder.lua
