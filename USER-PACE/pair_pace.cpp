/*
Copyright 2021 Yury Lysogorskiy^1, Cas van der Oord^2, Anton Bochkarev^1,
 Sarath Menon^1, Matteo Rinaldi^1, Thomas Hammerschmidt^1, Matous Mrovec^1,
 Aidan Thompson^3, Gabor Csanyi^2, Christoph Ortner^4, Ralf Drautz^1

^1: Ruhr-University Bochum, Bochum, Germany
^2: University of Cambridge, Cambridge, United Kingdom
^3: Sandia National Laboratories, Albuquerque, New Mexico, USA
^4: University of British Columbia, Vancouver, BC, Canada


    This FILENAME is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


//
// Created by Lysogorskiy Yury on 27.02.20.
//

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "pair_pace.h"
#include "atom.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "force.h"
#include "comm.h"
#include "memory.h"
#include "error.h"


#include "math_const.h"

#include "ace_version.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define MAXLINE 1024
#define DELTA 4

//added YL

//keywords for ACE evaluator style
#define RECURSIVE_KEYWORD "recursive"
#define PRODUCT_KEYWORD "product"


int elements_num_pace = 104;
char const *const elements_pace[104] = {"X", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na",
                                        "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn",
                                        "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr",
                                        "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb",
                                        "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd",
                                        "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir",
                                        "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
                                        "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr"
};

int AtomicNumberByName_pace(char *elname) {
    for (int i = 1; i < elements_num_pace; i++)
        if (strcmp(elname, elements_pace[i]) == 0)
            return i;
    return -1;
}


/* ---------------------------------------------------------------------- */
PairPACE::PairPACE(LAMMPS *lmp) : Pair(lmp) {
    //single_enable = 0;
    restartinfo = 0;
    one_coeff = 1;
    manybody_flag = 1;

    nelements = 0;

    ace = NULL;
    potential_file_name = NULL;
    elements = NULL;
    map = NULL;
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairPACE::~PairPACE() {
    if (copymode) return;

    if (elements)
        for (int i = 0; i < nelements; i++) delete[] elements[i];
    delete[] elements;


    delete[] potential_file_name;

    delete basis_set;
    delete ace;

    if (allocated) {
        memory->destroy(setflag);
        memory->destroy(cutsq);
        memory->destroy(map);
        memory->destroy(scale);
    }
}

/* ---------------------------------------------------------------------- */

void PairPACE::compute(int eflag, int vflag) {
    int i, j, ii, jj, inum, jnum;
    double delx, dely, delz, evdwl;
    double fij[3];
    int *ilist, *jlist, *numneigh, **firstneigh;

    ev_init(eflag, vflag);

    // downwards modified by YL

    double **x = atom->x;
    double **f = atom->f;
    tagint *tag = atom->tag;
    int *type = atom->type;

    // number of atoms in cell
    int nlocal = atom->nlocal;

    int newton_pair = force->newton_pair;

    // number of atoms including ghost atoms
    int nall = nlocal + atom->nghost;

    // inum: length of the neighborlists list
    inum = list->inum;

    // ilist: list of "i" atoms for which neighbor lists exist
    ilist = list->ilist;

    //numneigh: the length of each these neigbor list
    numneigh = list->numneigh;

    // the pointer to the list of neighbors of "i"
    firstneigh = list->firstneigh;

    if (inum != nlocal) {
        char str[128];
        snprintf(str,128,"inum: %d nlocal: %d are different",inum, nlocal);
        error->all(FLERR,str);
    }


    // Aidan Thompson told RD (26 July 2019) that practically always holds:
    // inum = nlocal
    // i = ilist(ii) < inum
    // j = jlist(jj) < nall
    // neighborlist contains neighbor atoms plus skin atoms,
    //       skin atoms can be removed by setting skin to zero but here
    //       they are disregarded anyway


    //determine the maximum number of neighbours
    int max_jnum = -1;
    int nei = 0;
    for (ii = 0; ii < list->inum; ii++) {
        i = ilist[ii];
        jnum = numneigh[i];
        nei = nei + jnum;
        if (jnum > max_jnum)
            max_jnum = jnum;
    }

    ace->resize_neighbours_cache(max_jnum);

    //loop over atoms
    for (ii = 0; ii < list->inum; ii++) {
        i = list->ilist[ii];
        const int itype = type[i];

        const double xtmp = x[i][0];
        const double ytmp = x[i][1];
        const double ztmp = x[i][2];

        jlist = firstneigh[i];
        jnum = numneigh[i];

        // checking if neighbours are actually within cutoff range is done inside compute_atom
        // mapping from LAMMPS atom types ('type' array) to ACE species is done inside compute_atom
        //      by using 'ace->element_type_mapping' array
        // x: [r0 ,r1, r2, ..., r100]
        // i = 0 ,1
        // jnum(0) = 50
        // jlist(neigh ind of 0-atom) = [1,2,10,7,99,25, .. 50 element in total]
        try {
            ace->compute_atom(i, x, type, jnum, jlist);
        } catch (exception &e) {
            error->all(FLERR, e.what());
            exit(EXIT_FAILURE);
        }
        // 'compute_atom' will update the `ace->e_atom` and `ace->neighbours_forces(jj, alpha)` arrays

        for (jj = 0; jj < jnum; jj++) {
            j = jlist[jj];
            const int jtype = type[j];
            j &= NEIGHMASK;
            delx = x[j][0] - xtmp;
            dely = x[j][1] - ytmp;
            delz = x[j][2] - ztmp;

            fij[0] = scale[itype][jtype]*ace->neighbours_forces(jj, 0);
            fij[1] = scale[itype][jtype]*ace->neighbours_forces(jj, 1);
            fij[2] = scale[itype][jtype]*ace->neighbours_forces(jj, 2);


            f[i][0] += fij[0];
            f[i][1] += fij[1];
            f[i][2] += fij[2];
            f[j][0] -= fij[0];
            f[j][1] -= fij[1];
            f[j][2] -= fij[2];

            // tally per-atom virial contribution
            if (vflag)
                ev_tally_xyz(i, j, nlocal, newton_pair, 0.0, 0.0,
                             fij[0], fij[1], fij[2],
                             -delx, -dely, -delz);
        }

        // tally energy contribution
        if (eflag) {
            // evdwl = energy of atom I
            evdwl = scale[1][1]*ace->e_atom;
            ev_tally_full(i, 2.0 * evdwl, 0.0, 0.0, 0.0, 0.0, 0.0);
        }
    }

    if (vflag_fdotr) virial_fdotr_compute();


    // end modifications YL
}

/* ---------------------------------------------------------------------- */

void PairPACE::allocate() {
    allocated = 1;
    int n = atom->ntypes;

    memory->create(setflag, n + 1, n + 1, "pair:setflag");
    memory->create(cutsq, n + 1, n + 1, "pair:cutsq");
    memory->create(map, n + 1, "pair:map");
    memory->create(scale, n + 1, n + 1,"pair:scale");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairPACE::settings(int narg, char **arg) {
    if (narg > 1) {
        error->all(FLERR,
                   "Illegal pair_style command. Correct form:\n\tpair_style pace\nor\n\tpair_style pace ");
        error->all(FLERR, RECURSIVE_KEYWORD);
        error->all(FLERR, "or\n\tpair_style pace ");
        error->all(FLERR, PRODUCT_KEYWORD);
    }
    recursive = true; // default evaluator style: RECURSIVE
    if (narg > 0) {
        if (strcmp(arg[0], RECURSIVE_KEYWORD) == 0)
            recursive = true;
        else if (strcmp(arg[0], PRODUCT_KEYWORD) == 0) {
            recursive = false;
        } else {
            error->all(FLERR,
                       "Illegal pair_style command: pair_style pace ");
            error->all(FLERR, arg[0]);
            error->all(FLERR, "\nCorrect form:\n\tpair_style pace\nor\n\tpair_style pace recursive");
        }
    }

    if (comm->me == 0) {
        if (screen) fprintf(screen, "ACE version: %d.%d.%d\n", VERSION_YEAR, VERSION_MONTH, VERSION_DAY);
        if (logfile) fprintf(logfile, "ACE version: %d.%d.%d\n", VERSION_YEAR, VERSION_MONTH, VERSION_DAY);

        if (recursive) {
            if (screen) fprintf(screen, "Recursive evaluator is used\n");
            if (logfile) fprintf(logfile, "Recursive evaluator is used\n");
        } else {
            if (screen) fprintf(screen, "Product evaluator is used\n");
            if (logfile) fprintf(logfile, "Product evaluator is used\n");
        }
    }


}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairPACE::coeff(int narg, char **arg) {

    if (narg < 4)
        error->all(FLERR,
                   "Incorrect args for pair coefficients. Correct form:\npair_coeff * * <potential.ace> elem1 elem2 ...");

    if (!allocated) allocate();

    //number of provided elements in pair_coeff line
    int ntypes_coeff = narg - 3;

    if (ntypes_coeff != atom->ntypes) {
        char error_message[1024];
        snprintf(error_message, 1024,
                 "Incorrect args for pair coefficients. You provided %d elements in pair_coeff, but structure has %d atom types",
                 ntypes_coeff, atom->ntypes);
        error->all(FLERR, error_message);
    }

    char *type1 = arg[0];
    char *type2 = arg[1];
    char *potential_file_name = arg[2];
    char **elemtypes = &arg[3];

    // insure I,J args are * *

    if (strcmp(type1, "*") != 0 || strcmp(type2, "*") != 0)
        error->all(FLERR, "Incorrect args for pair coefficients");


    //load potential file
    basis_set = new ACECTildeBasisSet();
    if (comm->me == 0) {
        if (screen) fprintf(screen, "Loading %s\n", potential_file_name);
        if (logfile) fprintf(logfile, "Loading %s\n", potential_file_name);
    }
    basis_set->load(potential_file_name);

    if (comm->me == 0) {
        if (screen) fprintf(screen, "Total number of basis functions\n");
        if (logfile) fprintf(logfile, "Total number of basis functions\n");

        for (SPECIES_TYPE mu = 0; mu < basis_set->nelements; mu++) {
            int n_r1 = basis_set->total_basis_size_rank1[mu];
            int n = basis_set->total_basis_size[mu];
            if (screen) fprintf(screen, "\t%s: %d (r=1) %d (r>1)\n", basis_set->elements_name[mu].c_str(), n_r1, n);
            if (logfile) fprintf(logfile, "\t%s: %d (r=1) %d (r>1)\n", basis_set->elements_name[mu].c_str(), n_r1, n);
        }
    }

    // read args that map atom types to pACE elements
    // map[i] = which element the Ith atom type is, -1 if not mapped
    // map[0] is not used

    ace = new ACERecursiveEvaluator();
    ace->set_recursive(recursive);
    ace->element_type_mapping.init(atom->ntypes + 1);

    for (int i = 1; i <= atom->ntypes; i++) {
        char *elemname = elemtypes[i - 1];
        int atomic_number = AtomicNumberByName_pace(elemname);
        if (atomic_number == -1) {
            char error_msg[1024];
            snprintf(error_msg, 1024, "String '%s' is not a valid element\n", elemname);
            error->all(FLERR, error_msg);
        }
        SPECIES_TYPE mu = basis_set->get_species_index_by_name(elemname);
        if (mu != -1) {
            if (comm->me == 0) {
                if (screen)
                    fprintf(screen, "Mapping LAMMPS atom type #%d(%s) -> ACE species type #%d\n", i, elemname, mu);
                if (logfile)
                    fprintf(logfile, "Mapping LAMMPS atom type #%d(%s) -> ACE species type #%d\n", i, elemname, mu);
            }
            map[i] = mu;
            ace->element_type_mapping(i) = mu; // set up LAMMPS atom type to ACE species  mapping for ace evaluator
        } else {
            char error_msg[1024];
            snprintf(error_msg, 1024, "Element %s is not supported by ACE-potential from file %s", elemname,
                     potential_file_name);
            error->all(FLERR, error_msg);
        }
    }

    // clear setflag since coeff() called once with I,J = * *
    int n = atom->ntypes;
    for (int i = 1; i <= n; i++) {
        for (int j = i; j <= n; j++) {
            setflag[i][j] = 1;
            scale[i][j] = 1.0;
        }
    }

    // set setflag i,j for type pairs where both are mapped to elements

    int count = 1;
    for (int i = 1; i <= n; i++)
        for (int j = i; j <= n; j++)
            if (map[i] >= 0 && map[j] >= 0) {
                setflag[i][j] = 1;
                count++;
            }

    if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients");

    ace->set_basis(*basis_set, 1);
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairPACE::init_style() {
    if (atom->tag_enable == 0)
        error->all(FLERR, "Pair style pACE requires atom IDs");
    if (force->newton_pair == 0)
        error->all(FLERR, "Pair style pACE requires newton pair on");

    // request a full neighbor list
    int irequest = neighbor->request(this, instance_me);
    neighbor->requests[irequest]->half = 0;
    neighbor->requests[irequest]->full = 1;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairPACE::init_one(int i, int j) {
    if (setflag[i][j] == 0) error->all(FLERR, "All pair coeffs are not set");
    //cutoff from the basis set's radial functions settings
    scale[j][i] = scale[i][j];
    return basis_set->radial_functions->cut(map[i], map[j]);
}

/* ---------------------------------------------------------------------- 
    extract method for extracting value of scale variable
 ---------------------------------------------------------------------- */
void *PairPACE::extract(const char *str, int &dim)
{
    dim = 2;
    if (strcmp(str,"scale") == 0) return (void *) scale;
    return NULL;
}

