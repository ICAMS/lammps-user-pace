# lammps-pyace Python Package

Minimal Python bindings for ML-PACE (Atomic Cluster Expansion) C++ library, focused on multispecies basis generation.

## Installation

### From source

```bash
cd /Users/mitch/github/lammps-pyace/python
pip install -e .
```

This will:
1. Build the ML-PACE C++ library
2. Compile the pybind11 Python extension `_basis`
3. Install the `lammps_pyace` Python package

### Requirements

- Python >= 3.8
- CMake >= 3.15
- C++17 compiler
- numpy >= 1.19.0
- pybind11 (automatically fetched if not found)

## Usage

```python
from lammps_pyace import (
    ACEBBasisSet,
    ACECTildeBasisSet,
    BBasisConfiguration,
    create_multispecies_basis_config
)

# Create a multispecies basis configuration
potential_config = {
    'elements': ['Al', 'Ni'],
    'deltaSplineBins': 0.001,
    'embeddings': {
        'ALL': {
            'ndensity': 1,
            'npot': 'FinnisSinclair',
            'fs_parameters': [1, 1],
            'rho_core_cut': 200000,
            'drho_core_cut': 250
        }
    },
    'bonds': {
        'ALL': {
            'rcut': 5.0,
            'dcut': 0.01,
            'radbase': 'ChebPow',
            'radparameters': [2.0],
            'NameOfCutoffFunction': 'cos'
        }
    },
    'functions': {
        'ALL': {
            'nradmax_by_orders': [5, 2, 2],
            'lmax_by_orders': [0, 1, 1]
        }
    }
}

# Create basis configuration
bbasis_config = create_multispecies_basis_config(potential_config)

# Create basis set
bbasis = ACEBBasisSet(bbasis_config)

# Save to file
bbasis_config.save('potential.yaml')
```

## Package Structure

```
lammps_pyace/
├── __init__.py                      # Main package exports
├── const.py                         # Constants
├── basisextension.py                # Basis extension utilities  
├── multispecies_basisextension.py  # Multispecies basis generation
├── data/
│   └── mus_ns_uni_to_rawlsLS_np_rank.pckl  # Precomputed basis data
└── _basis.so                        # Compiled C++ extension
```

## Building on the Cluster

The package uses CMake to build the C++ extension. When building on the cluster:

```bash
# Configure CMake build
cd /Users/mitch/github/lammps-pyace/python
mkdir build && cd build
cmake ..

# Build
cmake --build . -j4

# Install
pip install -e ..
```

## C++ Components

The package binds to:
- `../ML-PACE/ace/` - ACE basis function code
- `../ML-PACE/ace-evaluator/` - Evaluation infrastructure
- `helpers/` - Python binding helper functions

## Notes

- Package name uses underscores: `lammps_pyace`
- Distribution name can use hyphens: `lammps-pyace`
- All imports use: `from lammps_pyace import ...`
