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
	double c1=6.85e-5; 
	double c2=-2.49e-9; // These parameters come from Pivano et al. by trivial fitting

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
                        stor::Ex = c1*pow(x,2) + c2*pow(x,4);

                        // this is the analytical derivative of the potential dE/dx
                        stor::dEx= (2*c1*x + 4*c2*pow(x,3));
               
	return 0;

	}// end of function 
		
	// function which calculate the Domain wall width
	double calculate_DW(double phi){
		stor::Dw_size=Pi*(sqrt(2*stor::A/(stor::muMs*stor::Ms* sin(phi*one_rad)*sin(phi*one_rad) + stor::muMs*stor::H_demag))); 
	}// end of function calculate_DW

	// In this routine we calculate the Zeeman field taking into account the frequency of the field
	double Zeeman(double Dt){
		stor::V=stor::V0*sin(stor::omega*Dt*one_rad);
	
	}
}//end of namespace


