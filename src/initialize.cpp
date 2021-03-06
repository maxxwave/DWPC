// This file is dedicated to initialize the program variables from a file allowing the user to set preffered parameters
// This file takes part from DWPC code
//
//
// (c) author: Razvan Ababei
// University of Sheffield
// date: 11.08.2019
//
//

#include <cmath>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <stdio.h>
#include <string.h>
#include "../hdr/storage.h"
#include "../hdr/euler_integrator.h"
#include "../hdr/input_map.h"

namespace stor{
	void initialize(){
        // Input is read and stored in a map class
        input_map_t inputs;
        if( inputs.read_file("input") != 0 ) {
            std::cerr << "Error reading input file, exiting." << std::endl;
            exit(-1);
        }

        //std::cout << "Stored values: " << std::endl;
        //inputs.print();

        stor::Ms = inputs.get<double>("Ms");
        stor::L = inputs.get<double>("L");
        stor::Ly = inputs.get<double>("Ly");
        stor::Lz = inputs.get<double>("Lz");
        stor::A = inputs.get<double>("Aex");
        stor::alpha = inputs.get<double>("alpha");
		stor::V0 = inputs.get<double>("H");
        stor::freq = inputs.get<double>("f");
        stor::omega=2*Pi*stor::freq;
        integrate::Dt = inputs.get<double>("Dt");
        integrate::totaltime = inputs.get<double>("totaltime");
        integrate::out_time = inputs.get<double>("out_time");
        integrate::scheme = inputs.get<std::string>("integrator");
        stor::program = inputs.get<std::string>("program");

        std::cout<<"The program has been initialized with following parameters:"<<std::endl;
        std::cout<<"Saturation, Ms = "<<stor::Ms<<" A/m"<<std::endl;
        std::cout<<"Length of the strip = "<<stor::L<<" m"<<std::endl;
        std::cout<<"Width = "<<stor::Ly<<" m"<<std::endl;
        std::cout<<"Thickness = "<<stor::Lz<<" m"<<std::endl;
        std::cout<<"Exchange stiffness, A = "<<stor::A<<" J/m"<<std::endl;
        std::cout<<"Gilbert damping, alpha = "<<stor::alpha<<""<<std::endl;
        std::cout<<"Amplitude of oscillating field, H0 = "<<stor::V0<<" A/m"<<std::endl;
        std::cout<<"Frequency of the field, omega = " <<stor::omega<<" Hz"<<std::endl;
        std::cout<<"Integration time step, Dt = "<<integrate::Dt <<" s"<<std::endl;
        std::cout<<"Initialization completed!"<<std::endl;
        std::cout<<"=====================================================================<"<<std::endl;

        // Initialise the parameters for the integration
        integrate::prefac1 =(-stor::alpha*stor::gamma)/((1+pow(stor::alpha,2))*2*stor::Ms*stor::Lz*stor::Ly);
        integrate::prefac2 =stor::mu0*stor::gamma*stor::H_demag/2.0; //(2.0+2.0*stor::alpha*stor::alpha);
        integrate::prefac3 =-stor::gamma/((1+ stor::alpha*stor::alpha)*2*stor::Ms*stor::Lz*stor::Ly);
        integrate::prefac4 =-(stor::gamma*stor::alpha*stor::mu0*stor::H_demag)/(2+2*stor::alpha*stor::alpha);
        integrate::zeeman_prefac1 = stor::gamma*stor::mu0*stor::alpha/(stor::alpha*stor::alpha+1.0);
        integrate::zeeman_prefac2 = stor::gamma*stor::mu0/(1.0 + stor::alpha*stor::alpha);
	}
}//end of namespace
