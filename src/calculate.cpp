// Author: Razvan Ababei
// University of Sheffield
// date:20th of August, 2019
//
//
// This file is dedicated to create the chain including the potential for strip
// This model includes two anti-notches where the potential follows a form given by:
//
//
#include <iostream>
#include <vector>
#include <cmath>
#include "../hdr/storage.h"
//#include "../hdr/calculate.h"

namespace calculate{
	// Defining some prefactors where we incorporate the constants in order to not be called each time in the loop
	double prefac1 = (-stor::alpha*stor::gamma)/((1+pow(stor::alpha,2))*2*stor::Ms*stor::Lz*stor::Ly);
	double prefac2 = stor::mu0*stor::gamma*stor::H_demag/2.0; //(2.0+2.0*stor::alpha*stor::alpha);
	double prefac3 = -stor::gamma/((1+ stor::alpha*stor::alpha)*2*stor::Ms*stor::Lz*stor::Ly);
	double prefac4 = -(stor::gamma*stor::alpha*stor::mu0*stor::H_demag)/(2+2*stor::alpha*stor::alpha);
	double zeeman_prefac1 = stor::gamma*stor::mu0*stor::alpha/(stor::alpha*stor::alpha+1.0);
    double zeeman_prefac2 = stor::gamma*stor::mu0/(1.0 + stor::alpha*stor::alpha);

    const double one_rad=Pi/180.0;
	// Calculating the number of cells for a given L and cell_size
	int a = int(stor::L/stor::cell_size);
	int N=2*a+1;

	// some parameters from origin fit
	// The expression is F(X)= A0 + A1*X**2 + A2*X**2 + A3*X**3 + A4*X**4 + ... + A8*X**8
	// These values correspond to Py antinotches
	/*double a0 = 2.395e-20;
	double a1 = 6.741e-15;
	double a2 = -1.29e-7;
	double a3 = -2.94244;
	double a4 = -2.567e8;
	double a5 = 2.91623e14;
	double a6 = 6.77e21;
	double a7 = -7.728e27;
	double a8 = 2.6771e35;*/
	//These values correspond to Ni antinotches
	double a0 = 2.21117e-21;
	double a1 = -3.8271e-15;
	double a2 = -1.2866e-6;
	double a3 =  0.61164;
	double a4 =  1.632e8;
	double a5 = 0.0;
	double a6 = 0.0;
	double a7 = 0.0;
	double a8 = 0.0;


	// function which calculate the potential energy depending on x
	double update_energy_antinotches(double x){

		// this is the potential for notches
		// need to find out who is a and b ???
		//stor::Ex = a0 + a1*x + a2*x*x + a3*pow(x,3) + a4*pow(x,4) + a5*pow(x,5) + a6*pow(x,6) + a7*pow(x,7)+ a8*pow(x,8);

		// this is the analytical derivative of the potential dE/dx
		stor::dEx=a1+ 2*a2*x + 3*a3*x*x + 4*a4*pow(x,3) + 5*a5*pow(x,4) + 6*a6*pow(x,5) + 7*a7*pow(x,6) + 8*a8*pow(x,7);

		return 0;
	}// end of function

	// in this function we will calculate the pinning energy for anti-notches
	/*double update_energy_notches(double x){
		stor::Ex=U0-U1*(exp(-(double x + stor::xl)*(double x + stor::xl)/(stor::L*stor::L) )
			      -(exp(-(double x + stor::xl)*(double x + stor::xr)/(stor::L*stor::L) )
		stor::dEx=-2*U1/(stor::L*stor::L)*(x+stor::xl)*exp(-(x+stor::xl)*(x+stor::xl)/(stor::L*stor::L)) +
			  +2*U1/(stor::L*stor::L)*(x+stor::xr)*exp(-(x+stor::xr)*(x+stor::xr)/(stor::L*stor::L))

	}// end of function
	*/
	double update_energy();

	// function which calculate the Domain wall width
	double calculate_DW(double phi){
		stor::Dw_size=Pi*(sqrt(2*stor::A/(stor::muMs*stor::Ms* sin(phi)*sin(phi) + stor::muMs*stor::H_demag))); // Pivano form of DW
		//stor::Dw_size=sqrt(2*stor::A/(stor::mu0*stor::Ms*stor::Ms*(stor::Ny*sin(phi)*sin(phi) + stor::Nz*cos(phi)*cos(phi)))); // Matt form
		return 0;
	}// end of function calculate_DW

	// In this routine we calculate the Zeeman field taking into account the frequency of the field
	double Zeeman(double time){

		stor::V = stor::V0*sin(stor::omega*time);

		//this equation can be used for benchmark1 program
		//stor::V=stor::V0*cos(stor::omega*time);
		return 0 ;
	}

	// we define two function for speed and angular speed
	double phi_t(double dEx, double phi_rk, double H){
		return prefac3*dEx+prefac4*sin(2*phi_rk)+ zeeman_prefac2*H;
	}
	double x_t(double DWs, double phi, double phi_t){
		return prefac2*sin(2*phi)*DWs + stor::alpha*DWs*phi_t;
	}
}//end of namespace


