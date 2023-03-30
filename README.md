# lammps-user-pace

## Installation:

You could get the supported version of LAMMPS from [GitHub repository](https://github.com/lammps/lammps)

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


### Build for lammps compute

* First download lammps branch with PACE compute and FitSNAP

```
git clone -b compute-pace git@github.com:jmgoff/lammps_compute_PACE.git
cd lammps_compute_PACE
mkdir build && cd build
```

```
cmake -D LAMMPS_EXCEPTIONS=on -D PKG_PYTHON=on -D BUILD_SHARED_LIBS=on -D CMAKE_BUILD_TYPE=Debug -D PKG_ML-IAP=on -D PKG_ML-PACE=on -D PKG_ML-SNAP=on -D BUILD_MPI=on -D BUILD_OMP=off  -D CMAKE_INSTALL_PREFIX=<$HOME>/.local -D PKG_MOLECULE=on ../cmake/
```
* Next, download this modified lammps-user-pace repo that contains extra arrays for breaking out descriptor contributions

```
git clone git@github.com:jmgoff/lammps-user-pace-1.git
cp lammps-user-pace-1/ML-PACE/ace-evaluator/ace_evaluator.* ./lammps-user-pace-v.2022.09.27/ML-PACE/ace-evaluator/
```

```
make -j
make install
```

* Now, set up paths
```
INSTALL_PATH=CMAKE_INSTALL_PREFIX/.local
export PYTHONPATH=$PYTHONPATH:$INSTALL_PATH/lib/python3.<version>/site-packages
export PYTHONPATH=$PYTHONPATH:$INSTALL_PATH/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$INSTALL_PATH/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/<path>/<to>/<lammps>/build
```
