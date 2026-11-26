// DWPC - Domain Wall Propagation Code
// Copyright (C) 2019-2026 Razvan-Vasile Ababei and Matt Ellis
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This program is free software: you can redistribute it and/or modify it under the terms of
// the GNU General Public License as published by the Free Software Foundation, either version 3
// of the License, or (at your option) any later version. It is distributed WITHOUT ANY WARRANTY;
// see the file LICENSE for details. If you use DWPC in published work, please cite it as
// described in CITATION.cff.

#ifndef __STORAGE_H__
#define __STORAGE_H__

#include <iostream>
#include <vector>
#include <string>

#include "../hdr/arrays.h"


const double Pi =  3.1415926535897932384626;
const double Pi_sqr =  Pi*Pi;

namespace stor{

    extern std::string program;
	extern double L;
	extern double Lz, Ly;
	extern double Ms, mu0;
	extern double muMs;
	// Temperature dependence of the saturation magnetisation
	extern bool use_Ms_T;
	extern double Ms0, Tc, Ms_exponent, Ms_rescaling_exponent;
	extern double Aex0, Aex_exponent;
	extern double pinning_exponent, Ms_pin_ref, pin_scale;
	extern double Ms_T(double T);
	extern void update_magnetic_parameters();
	extern void compute_demag_factors();
	// Time dependent temperature (Gaussian pulse)
	extern bool use_T_pulse;
	extern double T0, T_pulse_peak, T_pulse_t0, T_pulse_sigma;
	extern double T_of_t(double time);
	extern void set_temperature(double T);
	extern void update_temperature(double time);
	// Electron and lattice temperatures (equal to T_sim unless the TTM is used)
	extern double T_e, T_l;
	// Appends "\tT_e\tT_l\tMs" (at the given time) to an output line
	extern void write_temperature_columns(std::ostream &out, double time);
	extern const char *temperature_columns_header;
	extern double Nx, Ny, Nz;
	extern double alpha;
	extern double A;
	extern double xl,xr;
	extern double V, V0, freq, omega;
    extern double H_const;
	extern double cell_size;
	extern int Nwires;
	extern std::vector<double> V0_mdw;
	extern std::vector<double> x_coord;
	extern std::vector<double> phi_coord;
	extern double T_sim;
	extern double gamma;
	//extern std::vector<double> H_demag;
	extern double H_demag;
	extern double Delta;
	extern double x_dw;
	extern double phi_dw;
	extern double Ex, dEx;
	extern double vx, phi_dt;
	extern double Dw_size;
	extern double u, mu_B, e_el, beta, P, j_dens;
	extern double rij;
	extern std::vector<double> H_DW;
	extern std::vector<double> H_dd;
	extern std::vector<double> j_dens_dw;
	extern std::vector<double> u_dw;
	extern double my;
	extern void initialize();
	extern double A0, A1, A2, A3, A4, A5, A6, A7, A8;
    extern bool use_DW_coupling;

    extern bool use_edge;
    extern double H_edge_max;
    extern double edge_scale;
    extern array_t<2,double> H_edge;
    extern array_t<2,double> E_edge;
} //end of namespace
#endif
