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
	const double rad=Pi/180.0;
	// Calculating the number of cells for a given L and cell_size
	int a = int(stor::L/stor::cell_size);
	int N=2*a+1;

	//defining the potential parameters
	double c1=1e-8, c2=1e-10;

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
		if(stor::x_coord[j]<=stor::xl || stor::x_coord[j]>=stor::xr){
			// this is the potential for notches
			// need to find out who is a and b ???
			stor::E_x[j] = stor::V*(c1*pow(stor::x_coord[j],2) + c2*pow(stor::x_coord[j],4));

			// this is the analytical derivative of the potential dE/dx
			stor::dE_x[j]= stor::V*(2*stor::x_coord[j] + 4*pow(stor::x_coord[j],3));
		}
		else{
			// this is the potential in between the anti-notches
			stor::E_x[j] = stor::V - stor::V*(exp(((pow(stor::x_coord[j]+stor::xl,2))/(4*stor::L*stor::L)))-
				       exp(-((pow(stor::x_coord[j]+stor::xr,2))/(4*stor::L*stor::L))));
			stor::dE_x[j] = stor::V*((-2*(stor::x_coord[j] + stor::xl)/stor::L)*exp(-((pow(stor::x_coord[j]+stor::xl,2))/(stor::L*stor::L)))-
                                       (-2*(stor::x_coord[j] + stor::xr)/stor::L)*exp(-((pow(stor::x_coord[j]+stor::xr,2))/(stor::L*stor::L))));

		}
	//std::cout<<stor::E_x[j]<<std::endl;

	}// end of for loop
	}// end of function initialize
	
	// function which calculate the potential energy depending on x
	double update_energy(double x){
		
		if(x<=stor::xl || x>=stor::xr){
                        // this is the potential for notches
                        // need to find out who is a and b ???
                        stor::Ex = stor::V*(c1*pow(x,2) + c2*pow(x,4));

                        // this is the analytical derivative of the potential dE/dx
                        stor::dEx= stor::V*(2*x + 4*pow(x,3));
                }
                else{
                        // this is the potential in between the anti-notches
                        stor::Ex = stor::V - stor::V*(exp(((pow(x+stor::xl,2))/(4*stor::L*stor::L)))-
                                       exp(-((pow(x+stor::xr,2))/(4*stor::L*stor::L))));
                        stor::dEx = stor::V*((-2*(x + stor::xl)/stor::L)*exp(-((pow(x+stor::xl,2))/(stor::L*stor::L)))-
                                       (-2*(x + stor::xr)/stor::L)*exp(-((pow(x+stor::xr,2))/(stor::L*stor::L))));

                }
	return 0;

	}// end of function 
		
	// function which calculate the Domain wall width
	double calculate_DW(){
		for (int i=0; i<stor::Dw.size(); i++){
			stor::Dw[i]=Pi*(sqrt(2*stor::A/(stor::muMs*stor::muMs* sin(stor::phi_coord[i]*rad)*sin(stor::phi_coord[i]*rad) + stor::muMs*stor::H_demag) )); 
		}//end of loop
	}// end of function calculate_DW
}//end of namespace


