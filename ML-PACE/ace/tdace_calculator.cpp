//
// Created by Yury Lysogorskiy on 24.06.24
//
#include "tdace_calculator.h"

void TDACECalculator::compute(ACEAtomicEnvironment &atomic_environment, bool compute_b_grad, bool verbose,
                              bool compute_projections) {
    evaluator.init_timers();
    evaluator.total_time_calc_timer.start();

    int i, j, jj;
    double fx, fy, fz, dx, dy, dz;

    energy = 0;


    energies.resize(atomic_environment.n_atoms_real);
    energies.fill(0);
    forces.resize(atomic_environment.n_atoms_real, 3);// per-atom forces
    forces.fill(0);

    virial.fill(0);


    //loop over atoms
    //determine the maximum number of neighbours
    int max_jnum = 0;
    for (i = 0; i < atomic_environment.n_atoms_real; ++i)
        if (atomic_environment.num_neighbours[i] > max_jnum)
            max_jnum = atomic_environment.num_neighbours[i];

    evaluator.resize_neighbours_cache(max_jnum);

    for (i = 0; i < atomic_environment.n_atoms_real; ++i) {
        evaluator.compute_atom(i,
                               atomic_environment.x,
                               atomic_environment.species_type,
                               atomic_environment.num_neighbours[i],
                               atomic_environment.neighbour_list[i]);
        //this will also update the e_atom and neighbours_forces(jj, alpha) array

        //update global energies and forces accumulators
        energies(i) = evaluator.e_atom;

        energy += evaluator.e_atom;


        const DOUBLE_TYPE xtmp = atomic_environment.x[i][0];
        const DOUBLE_TYPE ytmp = atomic_environment.x[i][1];
        const DOUBLE_TYPE ztmp = atomic_environment.x[i][2];

        for (jj = 0; jj < atomic_environment.num_neighbours[i]; jj++) {
            j = atomic_environment.neighbour_list[i][jj];

            dx = atomic_environment.x[j][0] - xtmp;
            dy = atomic_environment.x[j][1] - ytmp;
            dz = atomic_environment.x[j][2] - ztmp;

            fx = evaluator.neighbours_forces(jj, 0);
            fy = evaluator.neighbours_forces(jj, 1);
            fz = evaluator.neighbours_forces(jj, 2);

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
        }
    } // loop over atoms (i_at)

    evaluator.total_time_calc_timer.stop();

#ifdef FINE_TIMING
    if (verbose) {
        printf("(Calculator)Total time: %ld microseconds\n", evaluator.total_time_calc_timer.as_microseconds());
        printf("(Evaluator) Total time/at:    %ld microseconds\n",
               evaluator.per_atom_calc_timer.as_microseconds() / atomic_environment.n_atoms_real);

        printf("setup/atom: %ld microseconds\n",
               evaluator.setup_timer.as_microseconds() / atomic_environment.n_atoms_real);

        printf("A_construction/atom: %ld microseconds\n",
               evaluator.A_construction_timer.as_microseconds() / atomic_environment.n_atoms_real);

        printf("Energy/atom: %ld microseconds\n",
               evaluator.energy_calc_timer.as_microseconds() / atomic_environment.n_atoms_real);
        printf("Weights and theta/atom: %ld microseconds\n",
               evaluator.weights_and_theta_timer.as_microseconds() / atomic_environment.n_atoms_real);
        printf("Forces/atom: %ld microseconds\n",
               evaluator.forces_calc_loop_timer.as_microseconds() / atomic_environment.n_atoms_real);
    }
#endif


}

