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

// Created by Christoph Ortner  on 03.06.2020
// Adapted to use splines by Chuck Witt (April 2023)

#include "ace-evaluator/acejl_radial.h"

#include <functional>
#include <cmath>
#include <string>

using namespace std;

void ACEjlRadialFunctions::read_yaml(YAML_PACE::Node node) {

    auto yaml_map_bond_specifications = node["bonds"].as<map<vector<int>, YAML_PACE::Node>>();

    // loop over element pairs
    for (const auto &p: yaml_map_bond_specifications) {

        SPECIES_TYPE mu_i = p.first[0];
        SPECIES_TYPE mu_j = p.first[1];

        // extract spline information for this pair 
        auto bond_yaml = p.second;
        cut(mu_i,mu_j) = bond_yaml["rcut"].as<DOUBLE_TYPE>();
        int nbins = bond_yaml["nbins"].as<int>();
        auto splinenodalvals = bond_yaml["splinenodalvals"].as<map<int,vector<DOUBLE_TYPE>>>();
        auto splinenodalderivs = bond_yaml["splinenodalderivs"].as<map<int,vector<DOUBLE_TYPE>>>();
    
        // set up spline interpolator (see SplineInterpolator::setupSplines)
        // note that for SplineInterpolators, ntot=nlut is defined such 
        // that if r/delta>=ntot then f(r/delta)=0
        SplineInterpolator& spl = splines(mu_i, mu_j);
        spl.cutoff = cut(mu_i,mu_j);
        spl.ntot = nbins;
        spl.deltaSplineBins = spl.cutoff / spl.ntot;
        spl.num_of_functions = splinenodalvals.size();
        spl.values.resize(spl.num_of_functions);
        spl.derivatives.resize(spl.num_of_functions);
        spl.second_derivatives.resize(spl.num_of_functions);
        spl.nlut = spl.ntot;
        spl.rscalelookup = (DOUBLE_TYPE) spl.nlut / spl.cutoff;
        spl.invrscalelookup = 1.0 / spl.rscalelookup;
        spl.lookupTable.init(spl.ntot+1, spl.num_of_functions, 4);
        for (int n=0; n<spl.nlut-1; ++n) {
            for (int func_id=0; func_id<spl.num_of_functions; func_id++) {
                const auto d = spl.deltaSplineBins;
                DOUBLE_TYPE f0 = splinenodalvals[func_id][n];
                DOUBLE_TYPE f1 = splinenodalvals[func_id][n+1];
                DOUBLE_TYPE f0d1 = splinenodalderivs[func_id][n]*d;
                DOUBLE_TYPE f1d1 = splinenodalderivs[func_id][n+1]*d;
                // store c0 + c1*r + c2*r^2 + c3*r^3 coefficients
                spl.lookupTable(n,func_id,0) = f0;
                spl.lookupTable(n,func_id,1) = f0d1;
                spl.lookupTable(n,func_id,2) = 3.0*(f1-f0)-f1d1-2.0*f0d1;
                spl.lookupTable(n,func_id,3) = -2.0*(f1-f0)+f1d1+f0d1;
            }
        }
    }
}

void ACEjlRadialFunctions::setuplookupRadspline() {
}

void ACEjlRadialFunctions::init(
        NS_TYPE nradb, LS_TYPE lmax, NS_TYPE nradial,
        DOUBLE_TYPE deltaSplineBins, SPECIES_TYPE nelements,
        vector<vector<string>> radbasename) {

    // these spline members are unique to ACEjlRadialFunctions
    splines.init(nelements, nelements, "splines");

    // ... and the rest are inherited from AbstractRadialBasis
    this->nelements = nelements;
    cut.init(nelements, nelements, "cut");
    cut.fill(0.0);  // zeroed here, updated later
    dcut.init(nelements, nelements, "dcut");
    dcut.fill(1.0);

    this->deltaSplineBins = deltaSplineBins;
    this->lmax = lmax;
    this->nradial = nradial;
    nradbase = nradb;

    radbasenameij = radbasename;

    gr.init(nradb, "gr");
    dgr.init(nradb, "dgr");
    d2gr.init(nradbase, "d2gr");

    fr.init(nradial, lmax + 1, "fr");
    dfr.init(nradial, lmax + 1, "dfr");
    d2fr.init(nradial, lmax + 1, "d2fr");

    cr = 0.0;
    dcr = 0.0;

    crad.init(nelements, nelements, nradial, (lmax + 1), nradb, "crad");
    crad.fill(0.0);

    cut_in.init(nelements, nelements, "cut_in");
    cut_in.fill(0.0);

    dcut_in.init(nelements, nelements, "dcut_in");
    dcut_in.fill(1e-5);

    lambda.init(nelements, nelements, "lambda");
    lambda.fill(1.0);

    prehc.init(nelements, nelements, "prehc");
    prehc.fill(0.0);
    lambdahc.init(nelements, nelements, "lambdahc");
    lambdahc.fill(1.0);
}

void ACEjlRadialFunctions::evaluate(
        DOUBLE_TYPE r, NS_TYPE nradbase_c, NS_TYPE nradial_c,
        SPECIES_TYPE mu_i, SPECIES_TYPE mu_j,
        bool calc_second_derivatives) {

    auto &spline = splines(mu_i, mu_j);
    spline.calcSplines(r, calc_second_derivatives);

    for (NS_TYPE n = 0; n < gr.get_dim(0); n++) {
        gr(n) = spline.values(n);
        dgr(n) = spline.derivatives(n);
        if (calc_second_derivatives)
            d2gr(n) = spline.second_derivatives(n);
    }

    for (NS_TYPE n = 0; n < fr.get_dim(0); n++) {
        for (LS_TYPE l = 0; l < fr.get_dim(1); l++) {
            fr(n,l) = spline.values(n);
            dfr(n,l) = spline.derivatives(n);
            if (calc_second_derivatives)
                d2fr(n,l) = spline.second_derivatives(n);
        }
    }
}
