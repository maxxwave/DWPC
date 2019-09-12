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


int main(){	
	// In this step we initialize the variables from the input file set by the user
	stor::initialize();

	// declare the output file
	std::ofstream outputfile;
	outputfile.open("output");

	for (int l=-200;l<200; l++){
		calculate::update_energy(l*1e-9);
		//std::cout<<stor::Ex<<"\t"<<stor::dEx<<"\t"<<l*1e-9<<std::endl;
	}
	double time=0.0;
	// perform some integrations
	for (long int i=1; i<200000; i++){
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
return 1;	
	//close the file
	outputfile.close();
}// end of program 
