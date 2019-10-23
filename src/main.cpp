// Author: Razvan Ababei
// University of Sheffield
// date 20th of August, 2019
//
//
// This file is dedicated to create the main part of the Domain-Wall Propagation Code (DWPC)
//
//
//

#include <iostream>
#include <cstdlib>
#include <fstream>

// include the headers
#include "../hdr/storage.hdr"
#include "../hdr/calculate.hdr"
#include "../hdr/euler_integrator.hdr"
#include "../hdr/program.hdr"

int main(){	
	// In this step we initialize the variables from the input file set by the user
	stor::initialize();

	// declare the output file
	std::ofstream outputfile;
	outputfile.open("output");
	
	/*for (int l=-200;l<200; l++){
		calculate::update_energy_antinotches(l*1e-9);
		std::cout<<stor::Ex<<"\t"<<stor::dEx<<"\t"<<l*1e-9<<std::endl;
	}*/
	// calculate the bifurcation diagram
	//programs::bifurcation();

	//perform some equilibration steps
	double time=0.0;
	int Nsteps = integrate::totaltime / integrate::Dt;
	int Nout = 1;

	std::cout << "Using integrator: " << integrate::scheme << std::endl;
	std::cout << "Runtime = " << integrate::totaltime << ", Nsteps = " << Nsteps << std::endl;
	
	outputfile << "#time(s)      X       phi      dx/dt      V" << std::endl;
	// perform some integrations
	for (long int i=0; i<Nsteps; i++){
		for(long int j=0; j<Nout; j++){
			time += integrate::Dt;
			calculate::Zeeman(time);
			if( integrate::scheme.compare("EULER") == 0)
				integrate::euler();
			else if( integrate::scheme.compare("RK4") == 0)
				integrate::runge_kutta();
			else {
				std::cerr << "ERROR: Integrator not identified!" << std::endl;
				exit(-1);
			}
		}
		/*
	outputfile<<"x_dw="<<"\t"<<stor::x_dw<<"\t"
		<<"phi="<<"\t"<<stor::phi_dw<<"\t"
		<<"vx="<<"\t"<<stor::vx<<"\t"
		<<"V=" <<"\t"<<stor::V<<"\t"
		<<"t="<<"\t"<<i*integrate::Dt<<std::endl;*/
	outputfile << time << "\t" << stor::x_dw << "\t" << stor::phi_dw << "\t" << stor::vx << "\t" << stor::V <<std::endl;
	


	}		
	//close the file
	outputfile.close();
	std::cout<<"Simulation ended successfully!"<<std::endl;
	return 1;
}// end of program 
