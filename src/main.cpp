// DWPC - Domain Wall Propagation Code
// Copyright (C) 2019-2026 Razvan-Vasile Ababei and Matt Ellis
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This program is free software: you can redistribute it and/or modify it under the terms of
// the GNU General Public License as published by the Free Software Foundation, either version 3
// of the License, or (at your option) any later version. It is distributed WITHOUT ANY WARRANTY;
// see the file LICENSE for details. If you use DWPC in published work, please cite it as
// described in CITATION.cff.

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

    std::cout<<"DWPC - Domain Wall Propagation Code"<<std::endl;
    std::cout<<"Developers: Razvan-Vasile Ababei and Matt Ellis"<<std::endl;
    std::cout<<"University of Sheffield (2019); UAIC Iasi (2026). Licence: GPL-3.0-or-later"<<std::endl;
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
	    std::cout<<"Simulation ended successfully!"<<std::endl;
        return 0;
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

	if( stor::program.compare("RC_transient") == 0) {
        	std::cout << "Running Transient Reservoir Computing (RC) program." << std::endl;
        	reservoir::run_transient();
        	return 0;
    }

    if( stor::program.compare("multi_dw") == 0) {
        	std::cout << "MultiDW program." << std::endl;
        	reservoir::multi_wires_dynamics();
        	return 0;
    }

    if (stor::program.compare("spin_current") == 0) {
		std::cout<<"Running spin current program" <<std::endl;
		programs::spin_current1();
	}
    if (stor::program.compare("Lyapunov") == 0){
    		std::cout<<"Running Largest Lyapunov exponent calculations" << std::endl;
		programs::Lyapunov_exponent();
		return 0;
	
    }

	std::cout<<"Simulation ended successfully!"<<std::endl;
	return 1;
}// end of program
