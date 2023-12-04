Seeder
======

Octree mesh generator for Tree based Elemental Meshes.
The actual sources are found in the seeder-source repository,
which is included here in the `sdr` subdirectory.

Additionally all the other parts required for compilation of
the mesh generator are included in this repository as submodules.

Use `git clone --recurse-submodules` when cloning this repository
to fetch the gathered subdirectories from the various repositories.

Prerequisite for building the mesh generator is an installed Python,
Fortran compiler and MPI library.
For compilation you need to point `FC` to the appropiate MPI
compiler wrapper. (Usually `export FC=mpif90`).

The mesh generator can then be built with

```
bin/waf configure build
```

To install it, run:

```
bin/waf install
```

Run `bin/waf --help` to see all options.

Documentation
-------------

See the [documentation](https://geb.inf.tu-dresden.de/doxy/seeder/index.html)
for more details.
