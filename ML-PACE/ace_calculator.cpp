//
// Created by Yury Lysogorskiy on 13.03.2020.
//

#include "ace_calculator.h"

void ACECalculator::compute(ACEAtomicEnvironment &atomic_environment, bool verbose) {
    if (evaluator == nullptr) {
        throw std::invalid_argument("Evaluator is not set");
    }
    evaluator->init_timers();
    evaluator->total_time_calc_timer.start();

    int i, j, jj;
    double fx, fy, fz, dx, dy, dz;

    energy = 0;


    energies.resize(atomic_environment.n_atoms_real);
    energies.fill(0);
    forces.resize(atomic_environment.n_atoms_real, 3);// per-atom forces
    forces.fill(0);

    virial.fill(0);


    //loop over atoms
#ifdef PRINT_MAIN_STEPS
    printf("=====LOOP OVER ATOMS=====\n");
#endif
    //determine the maximum number of neighbours
    int max_jnum = -1;
    for (i = 0; i < atomic_environment.n_atoms_real; ++i)
        if (atomic_environment.num_neighbours[i] > max_jnum)
            max_jnum = atomic_environment.num_neighbours[i];

    evaluator->resize_neighbours_cache(max_jnum);

#ifdef EXTRA_C_PROJECTIONS
    basis_peratom_projections_rank1.resize(atomic_environment.n_atoms_real);
    basis_peratom_projections.resize(atomic_environment.n_atoms_real);
    gamma_grade.resize(atomic_environment.n_atoms_real);
#endif

    for (i = 0; i < atomic_environment.n_atoms_real; ++i) {
#ifdef PRINT_MAIN_STEPS
        printf("current central atom: %d\n", i);
#endif
        evaluator->compute_atom(i,
                                atomic_environment.x,
                                atomic_environment.species_type,
                                atomic_environment.num_neighbours[i],
                                atomic_environment.neighbour_list[i]);
        //this will also update the e_atom and neighbours_forces(jj, alpha) array
#ifdef EXTRA_C_PROJECTIONS
        basis_peratom_projections_rank1[i] = evaluator->basis_projections_rank1.to_vector();
        basis_peratom_projections[i] = evaluator->basis_projections.to_vector();
        gamma_grade[i] = evaluator->max_gamma_grade;
#endif

#ifdef DEBUG_FORCES_CALCULATIONS
        for (jj = 0; jj < atomic_environment.num_neighbours[i]; jj++) {
            printf("neighbour_forces(i=%d->j=%d)=(%f,%f,%f)\n", i,
                   atomic_environment.neighbour_list[i][jj],
                   evaluator->neighbours_forces(jj, 0),
                   evaluator->neighbours_forces(jj, 1),
                   evaluator->neighbours_forces(jj, 2)
                    );
        }
#endif
        //update global energies and forces accumulators
        energies(i) = evaluator->e_atom;

        energy += evaluator->e_atom;


        const DOUBLE_TYPE xtmp = atomic_environment.x[i][0];
        const DOUBLE_TYPE ytmp = atomic_environment.x[i][1];
        const DOUBLE_TYPE ztmp = atomic_environment.x[i][2];

        for (jj = 0; jj < atomic_environment.num_neighbours[i]; jj++) {
            j = atomic_environment.neighbour_list[i][jj];

            dx = atomic_environment.x[j][0] - xtmp;
            dy = atomic_environment.x[j][1] - ytmp;
            dz = atomic_environment.x[j][2] - ztmp;

            fx = evaluator->neighbours_forces(jj, 0);
            fy = evaluator->neighbours_forces(jj, 1);
            fz = evaluator->neighbours_forces(jj, 2);

            forces(i, 0) += fx;
            forces(i, 1) += fy;
            forces(i, 2) += fz;

            //virial f_dot_r, identical to LAMMPS virial_fdotr_compute
            virial(0) += dx * fx;
            virial(1) += dy * fy;
            virial(2) += dz * fz;
            virial(3) += dx * fy;
            virial(4) += dx * fz;
            virial(5) += dy * fz;

            // update forces only for real atoms
            if (j < atomic_environment.n_atoms_real) {
                forces(j, 0) -= fx;
                forces(j, 1) -= fy;
                forces(j, 2) -= fz;
            } else if (atomic_environment.origins != nullptr) { // map ghost j into true_j within periodic cell
                int true_j = atomic_environment.origins[j];
                if (true_j > atomic_environment.n_atoms_real)
                    throw invalid_argument(
                            "Inconsistency of atomic environment: origin index j = " + to_string(true_j) +
                            "out of real atom index range");
                forces(true_j, 0) -= fx;
                forces(true_j, 1) -= fy;
                forces(true_j, 2) -= fz;
            } else {
                throw invalid_argument(
                        "Atomic environment is not consistent: no origins array for mapping ghost atoms");
            }
#ifdef DEBUG_FORCES_CALCULATIONS
            printf("accumulated forces: F(i=%d)=(%f,%f,%f)\n", i, forces(i, 0), forces(i, 1), forces(i, 2));
            printf("accumulated forces: F(j=%d)=(%f,%f,%f)\n", j, forces(j, 0), forces(j, 1), forces(j, 2));
#endif
        }
    } // loop over atoms (i_at)

    evaluator->total_time_calc_timer.stop();

#ifdef FINE_TIMING
    if (verbose) {
        printf("   Total time: %ld microseconds\n", evaluator->total_time_calc_timer.as_microseconds());
        printf("Per atom time:    %ld microseconds\n",
               evaluator->per_atom_calc_timer.as_microseconds() / atomic_environment.n_atoms_real);


        printf("Loop_over_nei/atom: %ld microseconds\n",
               evaluator->loop_over_neighbour_timer.as_microseconds() / atomic_environment.n_atoms_real);

        printf("       Energy/atom: %ld microseconds\n",
               evaluator->energy_calc_timer.as_microseconds() / atomic_environment.n_atoms_real);

        printf("       Forces/atom: %ld microseconds\n",
               evaluator->forces_calc_loop_timer.as_microseconds() / atomic_environment.n_atoms_real);

//        printf("phi_recalcs/atom: %ld microseconds\n",
//               evaluator->phi_recalc_timer.as_microseconds() / atomic_environment.n_atoms_real);

        printf("     forces_neig: %ld microseconds\n",
               evaluator->forces_calc_neighbour_timer.as_microseconds() / atomic_environment.n_atoms_real);

    }
#endif


}

void ACECalculator::set_evaluator(ACEEvaluator &aceEvaluator) {
    this->evaluator = &aceEvaluator;
}
