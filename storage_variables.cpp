// Author: Razvan V Ababei
// University of Sheffield
// date:19th of August, 2019
//
//
// In this file we aim to declare and store the main variables of the code used throughout the modeling process
// This file is dedicated to store the variables for a chain whose magnetic properties come from a multiscale approach previously computed by a micromagnetic software
// The magnetic properties can be changed depending on the magnetic material. In this case we refer for Nickel.



#include <iostream>
#include <vector> 
#include "../hdr/storage.hdr"

namespace stor{

	// The length of the chain on x-component
	double  L=400.0;

	// Declare the section 
	double Lz=20.0, Ly=10.0;

	// Declare the saturation magnetisation measured in Tesla as B= mu_0 Ms
	double muMs=1.07;

	// Gilbert Damping
	double alpha=0.01;

	// Exchange Stiffness
	double A=1.3e-11; // J/m
	
	// Position of the notches or anti-notches. We define 2 nothces or anti-notches ditributted symmetrically from the middle of the chain
	// We define these as a ratio of the total length, L with respect to a symmetry 
	double xl=-0.75*L; 
	double xr=0.75*L;
	
	// The harmonic potential, V=V0cos(omega t)
	// In V0 we incorporate the Zeeman energy and the internal fields 
	// We note that Zeeman energy is given by an oscillatory field
	double V=10.0, V0=10.0;
	double omega=3.0; // in GHz

	//defining the length of the cell
	//needs to be around half of the exchange length
	double cell_size=3; // in nm

	// Defining the coordinates arrays
	// Within 1D model we assume two colective coordinates defined by the x and azimuthal angle
	std::vector <double> x_coord;
	std::vector <double> phi_coord;
	
	// defining the position of the DW
	double x_dw=0.0;
	double phi_dw=0.0;
	
	// defining the potential and the derivative 
	double Ex=0.0;
	double dEx=0.0;

	// Simulation temperature
	double T_sim=0.0;

	// Magnetic constants
	// Gyromagnetic factor
	double gamma=1.76e11; // in T^-1 s^-1
	// Magnetic permeability of free space
	double mu0=12.57e-7; // H/m

	// Demagnetization field arrays
	//std::vector <double> H_demag;
	double H_demag;

	// defining the domain width
	double Delta=0.0;
	
	// creating the potential arrays Ex and the derivative dEx	
	std::vector<double> E_x;
	std::vector<double> dE_x;
	
	// creating an array for the domain wall width
	std::vector<double> Dw;

}// end of namespace
