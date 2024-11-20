/*
 * Performant implementation of atomic cluster expansion and interface to LAMMPS
 *
 * Copyright 2021  (c) Yury Lysogorskiy^1, Cas van der Oord^2, Anton Bochkarev^1,
 * Sarath Menon^1, Matteo Rinaldi^1, Thomas Hammerschmidt^1, Matous Mrovec^1,
 * Aidan Thompson^3, Gabor Csanyi^2, Christoph Ortner^4, Ralf Drautz^1
 *
 * ^1: Ruhr-University Bochum, Bochum, Germany
 * ^2: University of Cambridge, Cambridge, United Kingdom
 * ^3: Sandia National Laboratories, Albuquerque, New Mexico, USA
 * ^4: University of British Columbia, Vancouver, BC, Canada
 *
 *
 * See the LICENSE file.
 * This FILENAME is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

// Created by Christoph Ortner on 03.06.2020
// Updated by Chuck Witt to use splines during April/May 2023

#ifndef ACEJL_RADIAL_FUNCTIONS_H
#define ACEJL_RADIAL_FUNCTIONS_H

#include "ace-evaluator/ace_arraynd.h"
#include "ace-evaluator/ace_types.h"
#include "ace-evaluator/ace_radial.h"
#include <yaml-cpp/yaml.h>

class ACEjlRadialFunctions : public AbstractRadialBasis {
public:

    Array2D<SplineInterpolator> splines;

//////////////////////////////////

    ACEjlRadialFunctions() = default;

    ~ACEjlRadialFunctions() override = default;

    void read_yaml(YAML_PACE::Node node);

    void init(NS_TYPE nradb, LS_TYPE lmax, NS_TYPE nradial,
              DOUBLE_TYPE deltaSplineBins, SPECIES_TYPE nelements,
              vector<vector<string>> radbasename) override;

    void evaluate(DOUBLE_TYPE r, NS_TYPE nradbase_c, NS_TYPE nradial_c,
                  SPECIES_TYPE mu_i, SPECIES_TYPE mu_j,
                  bool calc_second_derivatives = false) override;

    void setuplookupRadspline() override;

    ACEjlRadialFunctions *clone() const override {
        return new ACEjlRadialFunctions(*this);
    };

};

#endif
