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

        outputfile << "#time(ns)       X(nm)           phi             H0           H		v(m/s)" << std::endl;
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

                outputfile << std::setprecision(4)
                    << time*1e9 << "\t\t"
                    << std::setprecision(8)
                    << stor::x_dw*1e9 << "\t"
                    << stor::phi_dw << "\t"
                    << stor::V0 << "\t"
                    << stor::V << "\t"
                    << stor::vx <<std::endl;

            }
        }
        outputfile.close();
        return 1;

    }//end of benchmark1 function

    // in this subroutine we appply a oscillatory field and we extract the dynamics
    double benchmark2()

    {	// declare the output file
        std::ofstream outputfile;
        outputfile.open("output");
        //perform some equilibration steps
        double time=0.0;
        int Nsteps = std::round(integrate::totaltime / integrate::Dt);
        int Nout = std::round(integrate::out_time / integrate::Dt);

        std::cout << "Using integrator: " << integrate::scheme << std::endl;
        std::cout << "Runtime = " << integrate::totaltime << ", Nsteps = " << Nsteps << std::endl;
        std::cout << "Output time = " << integrate::out_time << ", Nout = " << Nout << std::endl;

        outputfile << "#time(ns)       X(nm)           phi             dx/dt           V" << std::endl;
        outputfile << std::fixed;

        double error = 0.0;

        // perform some integrations
        for (long int i=0; i<Nsteps/Nout; i++){
            for(long int j=0; j<Nout; j++){
                if( integrate::scheme.compare("EULER") == 0)
                    integrate::euler(time);
                else if( integrate::scheme.compare("RK4") == 0)
                    integrate::runge_kutta(time);
                else if( integrate::scheme.compare("RK45") == 0)
                    error = integrate::RK45::runge_kutta_45(time, integrate::Dt);
                else {
                    std::cerr << "ERROR: Integrator not identified!" << std::endl;
                    exit(-1);
                }
            }

            outputfile << std::fixed << std::setprecision(4)
                << time*1e9 << "\t\t"
                << std::setprecision(12)
                << stor::x_dw*1e9 << "\t"
                << stor::phi_dw << "\t"
                << stor::vx << "\t"
                << stor::V << "\t"
                << calculate::calculate_DW(stor::phi_dw)*1e9 << "\t"
                << std::scientific << error
                <<std::endl;

        }



    }// end of benchmark2

    double show_potential(){
        for (int l=-200;l<200; l++){
            calculate::update_energy_antinotches(l*1e-9);
            std::cout<<stor::Ex<<"\t"<<stor::dEx<<"\t"<<l*1e-9<<std::endl;
        }
    }//end of show_potential

    // in this subroutine we appply a oscillatory field and we extract the dynamics
    double benchmark3()

    {	// declare the output file
        std::ofstream outputfile;
        outputfile.open("Benchmark3.data");
        //perform some equilibration steps
        double time=0.0;
        int Nsteps = std::round(integrate::totaltime / integrate::Dt);
        int Nout = std::round(integrate::out_time / integrate::Dt);

        std::cout << "Using integrator: " << integrate::scheme << std::endl;
        std::cout << "Runtime = " << integrate::totaltime << ", Nsteps = " << Nsteps << std::endl;
        std::cout << "Output time = " << integrate::out_time << ", Nout = " << Nout << std::endl;

        outputfile << "#time(ns)       X(nm)           phi             dx/dt           V" << std::endl;
        outputfile << std::fixed;

        double error = 0.0;

        // perform some integrations
        for (long int i=0; i<Nsteps/Nout; i++){
            for(long int j=0; j<Nout; j++){
                if( integrate::scheme.compare("EULER") == 0)
                    integrate::euler(time);
                else if( integrate::scheme.compare("RK4") == 0)
                    integrate::multi_dw::runge_kutta( stor::x_coord, stor::phi_coord, time, integrate::Dt);
                else if( integrate::scheme.compare("RK45") == 0)
                    error = integrate::RK45::runge_kutta_45(time, integrate::Dt);
                else {
                    std::cerr << "ERROR: Integrator not identified!" << std::endl;
                    exit(-1);
                }
            }

            outputfile << std::fixed << std::setprecision(4)
                << time*1e9 << "\t\t"
                << std::setprecision(4)
                << stor::x_coord[0]*1e9 << "\t"
                << stor::phi_coord[0] << "\t"
                << stor::x_coord[1]*1e9 << "\t"
                << stor::phi_coord[1] << "\t"
                << stor::x_coord[2]*1e9 << "\t"
                << stor::phi_coord[2] << "\t"
                << std::scientific << error
                <<std::endl;

        }



    }// end of benchmark2

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

        outputfile << "#time(ns)       X(nm)           phi            J_dens           H		v(m/s)" << std::endl;
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

                outputfile << std::setprecision(4)
                    << time*1e9 << "\t\t"
                    << std::setprecision(8)
                    << stor::x_dw*1e9 << "\t"
                    << stor::phi_dw << "\t"
                    << stor::j_dens << "\t"
                    << stor::V << "\t"
                    << stor::vx <<std::endl;

            }
        }
        outputfile.close();
        return 1;
	
    }


}// end of namespace
