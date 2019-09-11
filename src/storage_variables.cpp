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
#include <cmath>
#include "../hdr/storage.hdr"
#define PI 3.14159265
namespace stor{

	// The length of the chain on x-component
	double  L=200.0e-9;

	// Declare the section 
	double Lz=50.0e-9, Ly=5.0e-9;
	
	//declare the magnetisation
	double Ms=860E3; // A/m
	double mu0=PI*4e-7; //T^2 J^-1 m^3

	// defining the Demag factors
	//double Nx=0.021829576;
	double Ny = 1 - (2/PI)*atan(L/Lz) + (1/(2*PI))*(Lz/L)*log10(1 + (L/Lz)*(L/Lz)) - (1/(2*PI))*(L/Lz)*log10(1 + (Lz/L)*(Lz/L));
	double Nz = 1 - (2/PI)*atan(Lz/L) + (1/(2*PI))*(L/Lz)*log10(1 + (Lz/L)*(Lz/L)) - (1/(2*PI))*(Lz/L)*log10(1 + (L/Lz)*(L/Lz));

	// Declare the saturation magnetisation measured in Tesla as B= mu_0 Ms
	double muMs=mu0*Ms;

	// Gilbert Damping
	double alpha=0.01;

	// Exchange Stiffness
	double A=1.3e-11; // J/m
	
	// The harmonic potential, V=V0cos(omega t)
	// In V0 we incorporate the Zeeman energy and the internal fields 
	// We note that Zeeman energy is given by an oscillatory field
	double V0=4000.0; // A/m which correspond to 18 Oe
	double V=0.0; // Instant field V=V0cos(wt)
	double omega=2*PI*400E6; // in Hz

	//defining the length of the cell
	//needs to be around half of the exchange length
	double cell_size=3e-10; // in m

	// Defining the coordinates arrays
	// Within 1D model we assume two colective coordinates defined by the x and azimuthal angle
	std::vector <double> x_coord;
	std::vector <double> phi_coord;
	
	// defining the position of the DW
	double x_dw=0.0;
	double phi_dw=1.57;
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
	double H_demag= Ms*(Nz-Ny);//Ms*(Nz-Ny); // This formula is taken from J. Dean et al., APL 2015

	// defining the domain width
	double Delta=0.0;
	
	// creating the potential arrays Ex and the derivative dEx	
	std::vector<double> E_x;
	std::vector<double> dE_x;
	
	// creating an array for the domain wall width
	std::vector<double> Dw;
	double Dw_size=0.0;

}// end of namespace
