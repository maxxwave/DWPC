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

    std::cout<<"Developers: Razvan Ababei, Matt Ellis and Tom Hayward"<<std::endl;
    std::cout<<"University of Sheffield"<<std::endl;
	// In this step we initialize the variables from the input file set by the user
    //

	stor::initialize();

	if(stor::program.compare("Potential")==0){
        std::cout<<"Running Potential program"<<std::endl;
		programs::show_potential();
	}

	if(stor::program.compare("Benchmark1")==0){
        std::cout<<"Running Benchmark-1 program"<<std::endl;
		programs::benchmark1();
	}
	if(stor::program.compare("Benchmark2")==0){
        std::cout<<"Running 'Benchmark-2 (Time-series) program"<<std::endl;
		programs::benchmark2();
	}
	if(stor::program.compare("Benchmark3")==0){
        std::cout<<"Running 'Benchmark-3 (Time-series) program"<<std::endl;
		programs::benchmark3();
	}
	if (stor::program.compare("Bifurcation")==0){
        std::cout<<"Running Bifurcation program"<<std::endl;
		programs::bifurcation();
	}

	if( stor::program.compare("RC") == 0) {
        	std::cout << "Running Reservoir Computing (RC) program." << std::endl;
        	reservoir::run();
        	return 0;
    }
	if (stor::program.compare("spoken_digit_recognition") == 0) {
		std::cout<<"Running spoken digits recognition program" <<std::endl;
		reservoir::run_spoken_recognition();
	}
	if( stor::program.compare("RC_field_sequence") == 0) {
        	std::cout << "Running Reservoir Computing (RC) program." << std::endl;
        	reservoir::run_field_sequence();
        	return 0;
    }


	std::cout<<"Simulation ended successfully!"<<std::endl;
	return 1;
}// end of program
