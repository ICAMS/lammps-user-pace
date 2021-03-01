# lammps-user-pace

## Installation:

### Build with `make`

1. Clone the repository using `https://github.com/ICAMS/lammps-user-pace.git` or download the zip files from [here](https://github.com/ICAMS/lammps-user-pace/archive/main.zip).
2. copy `USER-PACE` directory from the cloned repository into `lammps/src` folder.
3. Include `USER-PACE` in the compilation by running `make yes-user-pace` from the `lammps/src` folder.
4. Compile lammps as usual, i.e. `make serial`

### Build with `cmake`

1. Clone the repository using `https://github.com/ICAMS/lammps-user-pace.git` or download the zip files from [here](https://github.com/ICAMS/lammps-user-pace/archive/main.zip).
2. copy `USER-PACE` directory from the cloned repository into `lammps/src` folder.
3. Build LAMMPS using:
   ```
   cd lammps
   mkdir build
   cd build
   cmake PKG_USER-PACE=ON ../cmake
   make
   ```
   For more information see [here](https://lammps.sandia.gov/doc/Build_cmake.html).

### Install using `conda`

LAMMPS package can be installed using conda by `conda install -c conda-forge lammps`. This distribution includes the `USER-PACE` package.


