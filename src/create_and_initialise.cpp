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
#include "../hdr/storage.hdr"
//#include "../hdr/create.hdr"

namespace create{
	const double Pi=3.141592653589793238462643383279502884197169399375105820974944;
	const double one_rad=Pi/180.0;
	// Calculating the number of cells for a given L and cell_size
	int a = int(stor::L/stor::cell_size);
	int N=2*a+1;

	//defining the potential parameters
	double c1=-8.21e-7; 
	double c2=8.45e7; // These parameters come from Pivano et al. by trivial fitting
	
	// some parameters from origin fit
	// The expression is F(X)= A0 + A1*X**2 + A2*X**2 + A3*X**3 + A4*X**4 + ... + A8*X**8
	double a0 = 9.023e-20;	
	double a1 = -8.11e-14;
	double a2 = -1.30e-5;
	double a3 = 5.42;
	double a4 = 4.66e8;
	//double a5 = 1.068e15;	
	//double a6 = 1.62e22;
	//double a7 = -1.64e28;
	//double a8 = 1.32e35;

	// in this function we create the chain for a given L and cell_size		
	double create(){
	// Defining the number of cells
	// resizing arrays
	stor::x_coord.resize(N,0.0);
	stor::phi_coord.resize(N,0.0);
	stor::E_x.resize(N,0.0);
	stor::dE_x.resize(N,0.0);
	stor::Dw.resize(N,0.0);
	//stor::H_demag=1.0 //A/m
	}

	double initialize(){
	// initializing the coordinates arrays
	stor::x_coord[0]=-stor::L;
	for (int i=1; i<stor::x_coord.size(); i++){
		stor::x_coord[i] = stor::x_coord[i-1] + stor::cell_size;  
		// we set the angle to be zero ??
		stor::phi_coord[i]=0;
	}
	for (int j=0; j<stor::x_coord.size(); j++){
		// initialize the potential arrays
		// We distinguished two regime of the potential: out of the anti-notches and in between
		// this is the potential for notches
		// need to find out who is a and b ???
		stor::E_x[j] = (c1*pow(stor::x_coord[j],2) + c2*pow(stor::x_coord[j],4));

		// this is the analytical derivative of the potential dE/dx
		stor::dE_x[j]= (2*stor::x_coord[j] + 4*pow(stor::x_coord[j],3));
	
	//std::cout<<stor::E_x[j]<<std::endl;

	}// end of for loop
	}// end of function initialize
	
	// function which calculate the potential energy depending on x
	double update_energy(double x){
		
                        // this is the potential for notches
                        // need to find out who is a and b ???
                        stor::Ex =a0 + a1*x + a2*x*x + a3*pow(x,3) + a4*pow(x,4);// + a5*pow(x,5) + a6*pow(x,6) + a7*pow(x,7)+ a8*pow(x,8);

                        // this is the analytical derivative of the potential dE/dx
                        stor::dEx=a1+ 2*a2*x + 3*a3*x*x + 4*a4*pow(x,3);// + 5*a5*pow(x,4) + 6*a6*pow(x,5) + 7*a7*pow(x,6) + 8*a8*pow(x,7);

	}// end of function 
		
	// function which calculate the Domain wall width
	double calculate_DW(double phi){
		stor::Dw_size=Pi*(sqrt(2*stor::A/(stor::muMs*stor::Ms* sin(phi)*sin(phi) + stor::muMs*stor::H_demag))); 
	}// end of function calculate_DW

	// In this routine we calculate the Zeeman field taking into account the frequency of the field
	double Zeeman(double time){
		stor::V=stor::V0*sin(stor::omega*time);
	
	}
}//end of namespace


