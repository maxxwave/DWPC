// DWPC - Domain Wall Propagation Code
// Copyright (C) 2019-2026 Razvan-Vasile Ababei and Matt Ellis
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This program is free software: you can redistribute it and/or modify it under the terms of
// the GNU General Public License as published by the Free Software Foundation, either version 3
// of the License, or (at your option) any later version. It is distributed WITHOUT ANY WARRANTY;
// see the file LICENSE for details. If you use DWPC in published work, please cite it as
// described in CITATION.cff.

// Author: Razvan V Ababei
// University of Sheffield
// date:19.08.2019
// (c)
//
// In this file we aim to declare and store the main variables of the code used throughout the modeling process
// This file is dedicated to store the variables for a chain whose magnetic properties come from a multiscale approach previously computed by a micromagnetic software
// The magnetic properties can be changed depending on the magnetic material. In this case we refer for Py.



#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include "../hdr/storage.h"


namespace stor{

    std::string program("RC");
	double kb=1.38064852e-23;// m2Kgs-2K-1

	// The length of the chain on x-component
	double L=200.0e-9;

	// Declare the section
	double Ly=50e-9, Lz=5.0e-9;
	//Ly -width,w
	//Lz -thickness,t

	//declare the magnetisation
	double Ms = 477e3; // A/m
	double mu0 = Pi*4e-7; //T^2 J^-1 m^3

	// defining the Demag factors
	//double Nx=0.021829576;
	double Ny = 1 - (2/Pi)*atan(Ly/Lz) + (1/(2*Pi))*(Lz/Ly)*log(1 + (Ly/Lz)*(Ly/Lz)) - (1/(2*Pi))*(Ly/Lz)*log(1 + (Lz/Ly)*(Lz/Ly));
	double Nz = 1 - (2/Pi)*atan(Lz/Ly) + (1/(2*Pi))*(Ly/Lz)*log(1 + (Lz/Ly)*(Lz/Ly)) - (1/(2*Pi))*(Lz/Ly)*log(1 + (Ly/Lz)*(Ly/Lz));

	// Declare the saturation magnetisation measured in Tesla as B= mu_0 Ms
	double muMs=mu0*Ms;

	// Temperature dependence of Ms: Ms(T) = Ms0*(1 - (T/Tc)^Ms_rescaling_exponent)^Ms_exponent
	// Tc, Ms_exponent (beta) and Ms_rescaling_exponent (alpha) as fitted to atomistic
	// (VAMPIRE) M(T) curves, e.g. Evans et al., PRB 91, 144425 (2015), Eq. (5)
	bool use_Ms_T = false;
	double Ms0 = 477e3;      // A/m, saturation magnetisation at T=0 K
	double Tc = 0.0;         // K, Curie temperature
	double Ms_exponent = 0.0; // critical exponent of M(T)
	double Ms_rescaling_exponent = 1.0; // temperature rescaling exponent (1 = classical)

	// Exchange stiffness: A(T) = Aex0*(Ms(T)/Ms0)^Aex_exponent (Aex_exponent=0: constant)
	double Aex0 = 0.0;
	double Aex_exponent = 0.0;

	// Pinning energy E(x,T) = E(x)*(Ms(T)/Ms(T0))^pinning_exponent (0: no scaling)
	double pinning_exponent = 0.0;
	double Ms_pin_ref = 0.0;
	double pin_scale = 1.0;

	// Gaussian temperature pulse: T(t) = T0 + T_pulse_peak*exp(-(t-T_pulse_t0)^2/(2*T_pulse_sigma^2))
	bool use_T_pulse = false;
	double T0 = 0.0;            // K, base (bath) temperature
	double T_pulse_peak = 0.0;  // K, peak temperature rise above T0
	double T_pulse_t0 = 0.0;    // s, time of the pulse maximum
	double T_pulse_sigma = 0.0; // s, standard deviation of the Gaussian

	// Electron and lattice temperatures (two-temperature model)
	double T_e = 0.0;
	double T_l = 0.0;
	const char *temperature_columns_header = "\tT_e(K)\tT_l(K)\tMs(A/m)";

	// Gilbert Damping
	double alpha=0.02;

	// Exchange Stiffness
	double A=0.0; // J/m

	// The harmonic potential, V=V0cos(omega t)
	// In V0 we incorporate the Zeeman energy and the internal fields
	// We note that Zeeman energy is given by an oscillatory field
	double V0=0.0; // A/m which correspond to 18 Oe
	double V=0.0; // Instant field V=V0cos(wt)
	double freq=0.0;
	double omega=0.0;

    double H_const = 0.0; // constant field

	//defining the length of the cell
	//needs to be around half of the exchange length
	double cell_size=3e-10; // in m

	//declare the coefficients of the energy polynom
	double A0,A1,A2,A3,A4,A5,A6,A7,A8;

	//define the current parameters
	double e_el=1.602176634e-19;
	double beta=0.0;
	double j_dens=0.0;
	double P=0.0;
	double mu_B=9.274e-24;
	double u;

    //Number of wires
    int Nwires = 1;

    std::vector<double> V0_mdw;

	// Defining the coordinates arrays
	// Within 1D model we assume two colective coordinates defined by the x and azimuthal angle
	std::vector <double> x_coord;
	std::vector <double> phi_coord;

	// defining the position of the DW
	double x_dw=0.0;
	double phi_dw=0.0;

	double vx=0.0;
	double phi_dt=0.0;

	// defining the potential and the derivative
	double Ex=0.0;
	double dEx=0.0;

	// Simulation temperature
	double T_sim=0.0;

	// Magnetic constants
	// Gyromagnetic factor
	double gamma=1.76E11; // in T^-1 s^-1

	// Demagnetization field arrays
	//std::vector <double> H_demag;
	double H_demag= Ms*(Nz-Ny); // This formula is taken from J. Dean et al., APL 2015

	// defining the domain width
	double Delta=0.0;

	// creating the potential arrays Ex and the derivative dEx
	std::vector<double> E_x;
	std::vector<double> dE_x;
	// separation distance between nanowires
	double rij=50e-9;
	//std::vector<double> rij;
	std::vector<double> H_DW;
	std::vector<double> j_dens_dw;
	std::vector<double> u_dw;

	// creating an array for the domain wall width
	std::vector<double> Dw;
	std::vector<double> H_dd;
	std::vector<double> H_DD;
	double Dw_size=0.0;

	//declaring the my of the DW
	double my=7.1e-18; //SI

    bool use_DW_coupling = false;

    bool use_edge = false;

    double edge_scale=5.0e9;
    double H_edge_max = 0.0;
    array_t<2,double> H_edge;
    array_t<2,double> E_edge;


}// end of namespace
