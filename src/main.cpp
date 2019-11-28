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
#include <iomanip>
#include <cstdlib>
#include <fstream>
#include <cmath>

// include the headers
#include "../hdr/storage.h"
#include "../hdr/calculate.h"
#include "../hdr/euler_integrator.h"
#include "../hdr/program.h"
#include "../hdr/rc.h"

int main(){
	// In this step we initialize the variables from the input file set by the user
    //
    std::cout << integrate::prefac1 << std::endl;

	stor::initialize();
	std::cout<<stor::H_demag<<"\t"<<stor::Nz<<"\t"<<stor::Ny<<std::endl;
    std::cout << integrate::prefac1 << std::endl;

	if(stor::program.compare("Potential")==0){
		programs::show_potential();
	}

	if(stor::program.compare("Benchmark1")==0){
		programs::benchmark1();
	}
	if(stor::program.compare("Benchmark2")==0){
		programs::benchmark2();
	}
	if (stor::program.compare("Bifurcation")==0){
		programs::bifurcation();
	}

	if( stor::program.compare("RC") == 0) {
        	std::cout << "Running Reservoir Computing (RC) program." << std::endl;
        	reservoir::run();
        	return 0;
    	}


	std::cout<<"Simulation ended successfully!"<<std::endl;
	return 1;
}// end of program
