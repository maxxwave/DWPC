// This file is dedicated to initialize the program variables from a file allowing the user to set preffered parameters
// This file takes part from DWPC code
//
//
// (c) author: Razvan Ababei
// University of Sheffield
// date: 11.08.2019
//
//

#include <chrono>
#include <cmath>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <stdio.h>
#include <string.h>
#include "../hdr/storage.h"
#include "../hdr/calculate.h"
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
        stor::T_sim = inputs.get<double>("Temperature");
        stor::freq = inputs.get<double>("f");
        stor::P= inputs.get<double>("P");
        stor::j_dens= inputs.get<double>("j");
        stor::beta= inputs.get<double>("beta");
        stor::omega=2*Pi*stor::freq;
        integrate::Dt = inputs.get<double>("Dt");
        integrate::totaltime = inputs.get<double>("totaltime");
        integrate::out_time = inputs.get<double>("out_time");
        integrate::scheme = inputs.get<std::string>("integrator");
        stor::program = inputs.get<std::string>("program");

        stor::Nwires = inputs.get<int>("Nwires");
        stor::x_coord.assign(stor::Nwires, 0.0);
        stor::phi_coord.assign(stor::Nwires, 0.0);
        stor::V0_mdw.assign(stor::Nwires, stor::V0);
        stor::H_DW.assign(stor::Nwires, 0.0);

        //stor::x_coord[0] = -2.5e-8;
        //stor::x_coord[1] = 0.0;
        //stor::x_coord[2] = -2.5e-8;

        // initialize with polynomial coefficients (a0,a1, ..., a8)
        stor::A0=inputs.get<double>("a0", 0.0);
        stor::A1=inputs.get<double>("a1", 0.0);
        stor::A2=inputs.get<double>("a2", -1.29e-7);
        stor::A3=inputs.get<double>("a3", 0.0);
        stor::A4=inputs.get<double>("a4", 1.632e8);
        stor::A5=inputs.get<double>("a5", 0.0);
        stor::A6=inputs.get<double>("a6", 0.0);
        stor::A7=inputs.get<double>("a7", 0.0);
        stor::A8=inputs.get<double>("a8", 0.0);


        typedef std::chrono::high_resolution_clock myclock;
        //myclock::time_point beginning = myclock::now();
        myclock::time_point tp = myclock::now();
        // obtain a seed from the timer
        //myclock::duration d = myclock::now() - beginning;
        myclock::duration d = tp.time_since_epoch();
        calculate::random_seed = d.count();
        std::cout << calculate::random_seed << std::endl;
        calculate::seed_rng(calculate::random_seed);


        //if( stor::Nwires > 1 ){
        integrate::multi_dw::setup(stor::Nwires);
	    //stor::V0_mdw[1]=stor::V0;
	    //}

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
        std::cout<<"The coefficients of E(x) are:"<<"\t"<<stor::A0<<"\t"<<stor::A1<<"\t"<<stor::A2<<"\t"<<stor::A3<<"\t"<<stor::A4<<std::endl;
        std::cout<<"=====================================================================<"<<std::endl;

        // Initialise the parameters for the integration
        calculate::prefac1 =(-stor::alpha*stor::gamma)/((1+pow(stor::alpha,2))*2*stor::Ms*stor::Lz*stor::Ly);
        calculate::prefac2 =stor::mu0*stor::gamma*stor::H_demag/2.0; //(2.0+2.0*stor::alpha*stor::alpha);
        calculate::prefac3 =-stor::gamma/((1+ stor::alpha*stor::alpha)*2*stor::Ms*stor::Lz*stor::Ly);
        calculate::prefac4 =-(stor::gamma*stor::alpha*stor::mu0*stor::H_demag)/(2+2*stor::alpha*stor::alpha);
        calculate::zeeman_prefac1 = stor::gamma*stor::mu0*stor::alpha/(stor::alpha*stor::alpha+1.0);
        calculate::zeeman_prefac2 = stor::gamma*stor::mu0/(1.0 + stor::alpha*stor::alpha);
    }
}//end of namespace
