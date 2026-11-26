// DWPC - Domain Wall Propagation Code
// Copyright (C) 2019-2026 Razvan-Vasile Ababei and Matt Ellis
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This program is free software: you can redistribute it and/or modify it under the terms of
// the GNU General Public License as published by the Free Software Foundation, either version 3
// of the License, or (at your option) any later version. It is distributed WITHOUT ANY WARRANTY;
// see the file LICENSE for details. If you use DWPC in published work, please cite it as
// described in CITATION.cff.

#ifndef __CALCULATE_H__
#define __CALCULATE_H__

#include <vector>
#include "../hdr/storage.h"

namespace calculate{
    extern unsigned random_seed;
    void seed_rng(const unsigned);
	extern int N;
	extern int a;
	extern double rad;
	extern double c1,c2;
	extern double create();
	extern double initialize();
	extern double calculate_DW(double phi);
	extern double update_energy(double x);
	extern double update_energy_notches(double x);
	extern double update_energy_antinotches(double x);

    double compute_Hpin(double &x);

	extern double Zeeman(double Dt);
	extern double noise( double T, double DW );
	extern double current_prefac;
	extern double current(double time, double );
	extern double prefac1, prefac2, prefac3, prefac4, zeeman_prefac1, zeeman_prefac2;
	double x_t(double DWs, double phi, double phi_t);
	double phi_t(double dEx, double phi_rk, double H);
	extern double DW_coupling(std::vector <double> &);
    void gradient( double &, double &, double, double, const double);
    void gradient( std::vector<double> &, std::vector<double> &, std::vector<double> &, std::vector<double> &, const double);
    void BDF( double &, double);
    void BDF1( double &, double);
    void BDF2( double &, double);
    void BDF3( double &, double);
    void noise_gradient( double&, double&, double, double, double, double);
    void noise_gradient( std::vector<double>&, std::vector<double>&, std::vector<double> &, std::vector<double> &, std::vector<double> &, std::vector<double> &);
    double Uniform();
    double Normal();
    void Normal(std::vector<double> &);
}
#endif
