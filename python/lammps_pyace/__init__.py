"""
lammps-pyace: Minimal Python package for PACE multispecies basis generation

This package provides Python bindings to ML-PACE C++ library for creating
and manipulating ACE (Atomic Cluster Expansion) basis sets.
"""

# Import C++ extension classes
from lammps_pyace.basis import (
    # Core basis classes
    ACEBBasisSet,
    ACECTildeBasisSet,
    BBasisConfiguration,
    BBasisFunctionSpecification,
    BBasisFunctionsSpecificationBlock,
    
    # Radial functions
    ACERadialFunctions,
    
    # Helper functions
    Fexp,
    FexpShiftedScaled,
)

# Import Python utilities
from lammps_pyace.multispecies_basisextension import (
    create_multispecies_basis_config,
    extend_multispecies_basis,
    generate_species_keys,
)

from lammps_pyace.basisextension import (
    construct_bbasisconfiguration,
    extend_basis,
)

__all__ = [
    # Main classes
    'ACEBBasisSet',
    'ACECTildeBasisSet',
    'BBasisConfiguration',
    'BBasisFunctionSpecification',
    'BBasisFunctionsSpecificationBlock',
    'ACERadialFunctions',
    
    # Main functions
    'create_multispecies_basis_config',
    'construct_bbasisconfiguration',
    'extend_multispecies_basis',
    'extend_basis',
    'generate_species_keys',
    
    # Utility functions
    'Fexp',
    'FexpShiftedScaled',
]

__version__ = "0.1.0"
