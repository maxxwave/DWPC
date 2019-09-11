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
#include "../hdr/create.hdr"
#include "../hdr/euler_integrator.hdr"
//#include "euler_integrator.cpp"
//int main(int argc, char* argv[])

int main(){
	
	// declare the output file
	std::ofstream outputfile;
	outputfile.open ("output");
 
	// execute the creation function
	create::create();
	// execute the initialization function
	create::initialize();

	for (int l=-200;l<200; l++){
		create::update_energy(l*1e-9);
		std::cout<<stor::Ex<<"\t"<<stor::dEx<<"\t"<<l*1e-9<<std::endl;
	}
	double time=0.0;
	// perform some integrations
	for (long int i=1; i<200000; i++){
		for(long int j=1; j<10000; j++){
			time += integrate::Dt;
			create::Zeeman(time);
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
