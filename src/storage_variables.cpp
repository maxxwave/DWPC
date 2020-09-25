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
	double Ny = 1 - (2/Pi)*atan(Ly/Lz) + (1/(2*Pi))*(Lz/Ly)*log10(1 + (Ly/Lz)*(Ly/Lz)) - (1/(2*Pi))*(Ly/Lz)*log10(1 + (Lz/Ly)*(Lz/Ly));
	double Nz = 1 - (2/Pi)*atan(Lz/Ly) + (1/(2*Pi))*(Ly/Lz)*log10(1 + (Lz/Ly)*(Lz/Ly)) - (1/(2*Pi))*(Lz/Ly)*log10(1 + (Ly/Lz)*(Ly/Lz));

	// Declare the saturation magnetisation measured in Tesla as B= mu_0 Ms
	double muMs=mu0*Ms;

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

	//defining the length of the cell
	//needs to be around half of the exchange length
	double cell_size=3e-10; // in m


    //Number of wires
    int Nwires = 1;

    std::vector<double> V0_mdw;

	// Defining the coordinates arrays
	// Within 1D model we assume two colective coordinates defined by the x and azimuthal angle
	std::vector <double> x_coord;
	std::vector <double> phi_coord;

	// defining the position of the DW
	double x_dw=0;
	double phi_dw=0;

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

	// creating an array for the domain wall width
	std::vector<double> Dw;
	double Dw_size=0.0;

}// end of namespace
