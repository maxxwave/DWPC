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
	programs::bifurcation();

	//perform some equilibration steps
	double time=0.0;
	/*for (int k=0; k<stor::omega*integrate::Dt; k++){
			time += integrate::Dt;
			calculate::Zeeman(time);
			integrate::euler();
	}//end of equilibration loop
	*/
	// perform some integrations
	for (long int i=1; i<100000; i++){
		for(long int j=1; j<10000; j++){
			time += integrate::Dt;
			calculate::Zeeman(time);
			integrate::euler();
		}
	outputfile<<"x_dw="<<"\t"<<stor::x_dw<<"\t"
		<<"phi="<<"\t"<<stor::phi_dw<<"\t"
		<<"vx="<<"\t"<<stor::vx<<"\t"
		<<"V=" <<"\t"<<stor::V<<"\t"
		<<"t="<<"\t"<<10000*i*integrate::Dt<<std::endl;

	}		
	//close the file
	outputfile.close();
	return 1;
}// end of program 
