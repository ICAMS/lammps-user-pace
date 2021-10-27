# lammps-user-pace

## Installation:

Before the multispecies PACE will be merged into main branch of the official LAMMPS repository, you could get the unofficial version of LAMMPS from [here](https://github.com/yury-lysogorskiy/lammps)

### Build with `make`

Follow LAMMPS installation instructions

1. Go to `lammps/src` folder
2. Compile the ML-PACE library by running `make lib-pace args="-b"`
3. Include `ML-PACE` in the compilation by running `make yes-ml-pace`
4. Compile lammps as usual, i.e. `make serial` or `make mpi`.

### Build with `cmake`


1. Create build directory and go there with 

```
cd lammps
mkdir build
cd build
```

2. Configure the lammps build with

```
cmake -DCMAKE_BUILD_TYPE=Release -DPKG_ML-PACE=ON ../cmake 
```

or 

```
cmake -DCMAKE_BUILD_TYPE=Release -D BUILD_MPI=ON -DPKG_ML-PACE=ON ../cmake
```

For more information see [here](https://lammps.sandia.gov/doc/Build_cmake.html).

   
3. Build LAMMPS using `cmake --build .` or `make`

