//
// Created by Yury Lysogorskiy on 31.05.23.
//

#ifndef PYACE_TDACE_CALCULATOR_H
#define PYACE_TDACE_CALCULATOR_H

#include "tdace_evaluator.h"
#include "extra/ace_atoms.h"

class TDACECalculator {
    TDACEBEvaluator &evaluator;
public:
    //total energy of ACEAtomicEnvironment
    DOUBLE_TYPE energy = 0;
    //total forces array
    //forces(i,3), i = 0..num_of_atoms-1
    Array2D<DOUBLE_TYPE> forces = Array2D<DOUBLE_TYPE>("forces");

    //stresses
    Array1D<DOUBLE_TYPE> virial = Array1D<DOUBLE_TYPE>(6, "virial");

    //Per-atom energies
    //energies(i), i = 0..num_of_atoms-1
    Array1D<DOUBLE_TYPE> energies = Array1D<DOUBLE_TYPE>("energies");

    explicit TDACECalculator(TDACEBEvaluator &aceEvaluator) : evaluator(aceEvaluator) {};


    //compute the energies and forces for each atoms in atomic_environment
    //results are stored in forces and energies arrays
    void compute(ACEAtomicEnvironment &atomic_environment, bool compute_b_grad = false, bool verbose = false, bool compute_projections = false);

};

#endif //PYACE_TDACE_CALCULATOR_H
