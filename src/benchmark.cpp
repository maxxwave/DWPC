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
// date:  11.11.2019
// This file is dedicated to perform some benchmark programs of DWPC
//
//
//



#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include "../hdr/calculate.h"
#include "../hdr/storage.h"
#include "../hdr/euler_integrator.h"
#include "../hdr/program.h"



namespace programs{
    //                                                                                           ____
    // Benchmark 1: This function will solve the dynamics of the DW for a square input field ___|
    double benchmark1()
    {
        // declare the output file
        std::ofstream outputfile;
        outputfile.open("Benchmark.data");
        double time=0.0;
        int Nsteps = std::round(integrate::totaltime / integrate::Dt);
        int Nout = std::round(integrate::out_time / integrate::Dt);
        double no_h_values=2;

        outputfile << "#time(ns)       X(nm)           phi             H0           H		v(m/s)" << stor::temperature_columns_header << std::endl;
        outputfile << std::fixed;

        double error;
        for (int k=0; k<no_h_values;k++)
        {
            // update the field
            stor::V0 *= (k+1);

            // perform some integrations
            for (long int i=0; i<Nsteps/Nout; i++)
            {

                for(long int j=0; j<Nout; j++)
                {
                    if(integrate::scheme.compare("EULER") == 0)
                        integrate::euler(time);
                    else if( integrate::scheme.compare("RK4") == 0)
                        integrate::runge_kutta(time);
                   // else if( integrate::scheme.compare("RK45") == 0)
                   //     integrate::RK45::Runge_Kutta_45(time, integrate::Dt);
                    else {
                        std::cerr << "ERROR: Integrator not identified!" << std::endl;
                        exit(-1);
                    }
                }

                outputfile << std::setprecision(10)
                    << time*1e9 << "\t\t"
                    << std::setprecision(8)
                    << stor::x_dw*1e9 << "\t"
                    << stor::phi_dw << "\t"
                    << stor::V0 << "\t"
                    << stor::V << "\t"
                    << stor::vx;
                stor::write_temperature_columns(outputfile, time);
                outputfile << std::endl;

            }
        }
        outputfile.close();
        return 1;

    }//end of benchmark1 function

    // in this subroutine we appply a oscillatory field and we extract the dynamics
    double benchmark2()

    {
        // declare the output file
        std::ofstream outputfile;
        outputfile.open("output");
        //perform some equilibration steps
        double time=0.0;
        int Nsteps = std::round(integrate::totaltime / integrate::Dt);
        int Nout = std::round(integrate::out_time / integrate::Dt);

        std::cout << "Using integrator: " << integrate::scheme << std::endl;
        std::cout << "Runtime = " << integrate::totaltime << ", Nsteps = " << Nsteps << std::endl;
        std::cout << "Output time = " << integrate::out_time << ", Nout = " << Nout << std::endl;

        outputfile << "#time(ns)       X(nm)           phi             dx/dt (m/s)           V" << stor::temperature_columns_header << std::endl;
        outputfile << std::fixed;

        double error = 0.0;
        int loop_steps = std::round(Nsteps / Nout);

        // perform some integrations
        for (long int i=0; i<loop_steps; i++){
            for(long int j=0; j<Nout; j++){
                if( integrate::scheme.compare("EULER") == 0)
                    integrate::euler(time);
                else if( integrate::scheme.compare("RK4") == 0)
                    integrate::runge_kutta(time);
                else if( integrate::scheme.compare("HEUN") == 0)
                    integrate::heun(time);
                else if( integrate::scheme.compare("RK45") == 0)
                    error = integrate::RK45::runge_kutta_45(time, integrate::Dt);
                else {
                    std::cerr << "ERROR: Integrator not identified!" << std::endl;
                    exit(-1);
                }
            }

            outputfile << std::fixed << std::setprecision(10)
                << time*1e9 << "\t\t"
                << std::setprecision(12)
                << stor::x_dw*1e9 << "\t"
                << stor::phi_dw << "\t"
                << stor::vx << "\t"
                << stor::V << "\t"
                << calculate::current(time, stor::j_dens) << "\t"
                << calculate::calculate_DW(stor::phi_dw)*1e9 << "\t"
                << std::scientific << error;
            stor::write_temperature_columns(outputfile, time);
            outputfile << std::endl;

        }

        std::cerr << "Finished benchmark 2" <<std::endl;
        outputfile.close();
        return 1;

    }// end of benchmark2

    double show_potential(){
        for (int l=-200;l<200; l++){
            calculate::update_energy_antinotches(l*1e-9);
            std::cout<<stor::Ex<<"\t"<<stor::dEx<<"\t"<<l*1e-9<<std::endl;
        }
    }//end of show_potential

    // in this subroutine we appply a oscillatory field and we extract the dynamics
    // Time series of Nwires coupled DWs: x_i(t), phi_i(t) for every wire
    double benchmark3()
    {
        std::ofstream outputfile;
        outputfile.open("Benchmark3.data");
        double time=0.0;
        long int Nsteps = std::lround(integrate::totaltime / integrate::Dt);
        long int Nout = std::lround(integrate::out_time / integrate::Dt);

        if( integrate::scheme.compare("RK4") != 0 ){
            std::cerr << "ERROR: Benchmark3 (multi-DW) supports integrator=RK4 only." << std::endl;
            exit(-1);
        }

        std::cout << "Using integrator: " << integrate::scheme << std::endl;
        std::cout << "Runtime = " << integrate::totaltime << ", Nsteps = " << Nsteps << std::endl;
        std::cout << "Output time = " << integrate::out_time << ", Nout = " << Nout << std::endl;

        outputfile << "#time(ns)";
        for( int k = 0; k < stor::Nwires; k++)
            outputfile << "\tX" << k << "(nm)\tphi" << k;
        outputfile << "\tH(A/m)" << stor::temperature_columns_header << std::endl;

        // initial state
        long int loop_steps = Nsteps / Nout;
        for (long int i=0; i<=loop_steps; i++){
            if( i > 0 ){
                for(long int j=0; j<Nout; j++)
                    integrate::multi_dw::runge_kutta( stor::x_coord, stor::phi_coord, time, integrate::Dt);
            }
            outputfile << std::fixed << std::setprecision(10) << time*1e9;
            outputfile << std::setprecision(6);
            // initial state is written before the first step
            for( int k = 0; k < stor::Nwires; k++)
                outputfile << "\t" << stor::x_coord[k]*1e9 << "\t" << stor::phi_coord[k];
            outputfile << "\t" << stor::V0_mdw[0]*sin(stor::omega*time);
            stor::write_temperature_columns(outputfile, time);
            outputfile << std::endl;
        }
        outputfile.close();
        return 1;
    }// end of benchmark3

   // double multi_DW(){

    //}
    //

    // In this routine the dynamics of DW will be calculated while applying a square current pulse
    // _____|----|_____
    //
    double spin_current1(){
	// declare the output file
        std::ofstream outputfile;
        outputfile.open("output.data");
        double time=0.0;
        int Nsteps = std::round(integrate::totaltime / integrate::Dt);
        int Nout = std::round(integrate::out_time / integrate::Dt);
        double no_h_values=10;

        outputfile << "#time(ns)       X(nm)           phi            J_dens           H		v(m/s)" << stor::temperature_columns_header << std::endl;
        outputfile << std::fixed;

	    double const amp=stor::j_dens;
        double error;

        for (int k=0; k<no_h_values;k++)
        {
            stor::j_dens = pow(-1,k)*amp;

	    for (long int i=0; i<Nsteps/Nout; i++){

		for(long int j=0; j<Nout; j++)
                {
                    if(integrate::scheme.compare("EULER") == 0)
                        integrate::euler(time);
                    else if( integrate::scheme.compare("RK4") == 0)
                        integrate::runge_kutta(time);
                   // else if( integrate::scheme.compare("RK45") == 0)
                    //    integrate::RK45::Runge_Kutta_45(time, integrate::Dt);
                    else {
                        std::cerr << "ERROR: Integrator not identified!" << std::endl;
                        exit(-1);
                    }
                }

                outputfile << std::setprecision(10)
                    << time*1e9 << "\t\t"
                    << std::setprecision(8)
                    << stor::x_dw*1e9 << "\t"
                    << stor::phi_dw << "\t"
                    << stor::j_dens << "\t"
                    << stor::V << "\t"
                    << stor::vx;
                stor::write_temperature_columns(outputfile, time);
                outputfile << std::endl;

            }
        }
        outputfile.close();
        return 1;

    }

double Lyapunov_exponent(){
	//in this vector we store the displacements x0,x1, .., xn
	std::vector<double> data_series;
	//in this vector we store the lyapunov values
	std::vector<double> lyapunov_series;


	// declare the output file
        std::ofstream outputfile;
        outputfile.open("output");
        //perform some equilibration steps
        double time=0.0;
        int Nsteps = std::round(integrate::totaltime / integrate::Dt);
        int Nout = std::round(integrate::out_time / integrate::Dt);

        std::cout << "Using integrator: " << integrate::scheme << std::endl;
        std::cout << "Runtime = " << integrate::totaltime << ", Nsteps = " << Nsteps << std::endl;
        std::cout << "Output time = " << integrate::out_time << ", Nout = " << Nout << std::endl;

        outputfile << "#time(ns)       X(nm)           phi             dx/dt (m/s)           V" << stor::temperature_columns_header << std::endl;
        outputfile << std::fixed;

        double error = 0.0;
        int loop_steps = std::round(Nsteps / Nout);

        // perform some integrations
        for (long int i=0; i<loop_steps; i++){
            for(long int j=0; j<Nout; j++){
                if( integrate::scheme.compare("EULER") == 0)
                    integrate::euler(time);
                else if( integrate::scheme.compare("RK4") == 0)
                    integrate::runge_kutta(time);
                else if( integrate::scheme.compare("HEUN") == 0)
                    integrate::heun(time);
                else if( integrate::scheme.compare("RK45") == 0)
                    error = integrate::RK45::runge_kutta_45(time, integrate::Dt);
                else {
                    std::cerr << "ERROR: Integrator not identified!" << std::endl;
                    exit(-1);}//end of else loop
            }//end of Nout loop

            outputfile << std::fixed << std::setprecision(10)
                << time*1e9 << "\t\t"
                << std::setprecision(12)
                << stor::x_dw*1e9 << "\t"
                << stor::phi_dw << "\t"
                << stor::vx << "\t"
                << stor::V << "\t"
                << calculate::current(time, stor::j_dens) << "\t"
                << calculate::calculate_DW(stor::phi_dw)*1e9 << "\t"
                << std::scientific << error;
            stor::write_temperature_columns(outputfile, time);
            outputfile << std::endl;
	    //std::cout<<"value of field is ....."<<stor::V<<std::endl;
	    //if(stor::V0-0.0001<stor::V){
		//std::cout<<"value of V is...."<<stor::V<<std::endl;
	    	data_series.push_back(stor::x_dw); 
	    //}

        }// end of integration loop	

	std::vector<double> data_series2;
	// reset the initial position 
	stor::x_dw=1e-10;
	stor::phi_dw=1e-11;
	time=0.0;
	//data_series2[0]=5e-9;
	// perform some integrations
        for (long int i=0; i<loop_steps; i++){
            for(long int j=0; j<Nout; j++){
                if( integrate::scheme.compare("EULER") == 0)
                    integrate::euler(time);
                else if( integrate::scheme.compare("RK4") == 0)
                    integrate::runge_kutta(time);
                else if( integrate::scheme.compare("HEUN") == 0)
                    integrate::heun(time);
                else if( integrate::scheme.compare("RK45") == 0)
                    error = integrate::RK45::runge_kutta_45(time, integrate::Dt);
                else {
                    std::cerr << "ERROR: Integrator not identified!" << std::endl;
                    exit(-1);}//end of else loop
            }//end of Nout loop

            outputfile << std::fixed << std::setprecision(10)
                << time*1e9 << "\t\t"
                << std::setprecision(12)
                << stor::x_dw*1e9 << "\t"
                << stor::phi_dw << "\t"
                << stor::vx << "\t"
                << stor::V << "\t"
                << calculate::current(time, stor::j_dens) << "\t"
                << calculate::calculate_DW(stor::phi_dw)*1e9 << "\t"
                << std::scientific << error;
            stor::write_temperature_columns(outputfile, time);
            outputfile << std::endl;
	    	data_series2.push_back(stor::x_dw); 
	    //}

        }// end of integration loop

	std::cout<<"Length of data_series is ... "<<data_series.size()<<std::endl;
        std::cerr << "Dynamics calculations have been finished" <<std::endl;
        outputfile.close();

	std::cout<<"Calculating the Lyapunov series .. "<< std::endl;
	double lyapunov_sum=0.0;
	double delta_x0=std::max(fabs(data_series2[0] - data_series[0]),1e-11);

	for (int i=0; i<data_series.size()-1; i++){
		double delta_x = std::max(fabs(data_series2[i]-data_series[i]),1e-23); // /(integrate::Dt*Nout);
		//double delta_x = fabs(data_series2[i]-data_series[i]); // /(integrate::Dt*Nout);
		std::cout<<delta_x<<std::endl;
		lyapunov_series.push_back(log(fabs(delta_x/delta_x0))/Nout);
		//lyapunov_series.push_back(log(delta_x / delta_x0) / (integrate::Dt * Nout));	
	}
	
	//std::cout<<"incepe afisarea ...................."<<std::endl;
	for (int i=0; i<lyapunov_series.size(); i++){
		lyapunov_sum += lyapunov_series[i];
		//std::cout<<lyapunov_series.size()<<"\t"<<data_series.size()<<std::endl;
		std::cout<<std::fixed<<std::setprecision(4)<<
			lyapunov_sum<<"\t"<<
			lyapunov_sum/(i+1)<<"\t"<<
			lyapunov_series[i]<<"\t"
			<<i<<std::scientific<<std::endl;
	}
	std::cout<<"delta x0 is   "<<delta_x0<<"delta x2   "<<"Last Lyap"<< lyapunov_series[lyapunov_series.size()-1]<<std::endl;	
        return 1;

}

}// end of namespace
