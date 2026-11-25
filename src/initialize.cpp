// This file is dedicated to initialize the program variables from a file allowing the user to set preffered parameters
// This file takes part from DWPC code
//
//
// (c) author: Razvan Ababei
// University of Sheffield
// date: 11.08.2019
//
//

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <iomanip>
#include <stdio.h>
#include <string.h>
#include "../hdr/storage.h"
#include "../hdr/calculate.h"
#include "../hdr/euler_integrator.h"
#include "../hdr/input_map.h"
#include "../hdr/arrays.h"
#include "../hdr/ttm.h"

namespace stor{

    // Saturation magnetisation at temperature T (in K):
    //   Ms(T) = Ms0 * (1 - T/Tc)^Ms_exponent,   for T < Tc
    //   Ms(T) = 0,                              for T >= Tc
    // If use_Ms_T is false the constant Ms0 is returned.
    double Ms_T(double T){
        if( !stor::use_Ms_T ) return stor::Ms0;
        if( T <= 0.0 ) return stor::Ms0;
        if( T >= stor::Tc ) return 0.0;
        return stor::Ms0 * pow(1.0 - T/stor::Tc, stor::Ms_exponent);
    }

    // Demagnetising factors of the strip cross-section (depend only on Ly, Lz)
    void compute_demag_factors(){
        stor::Ny = 1 - (2/Pi)*atan(stor::Ly/stor::Lz) + (1/(2*Pi))*(stor::Lz/stor::Ly)*log10(1 + (stor::Ly/stor::Lz)*(stor::Ly/stor::Lz)) - (1/(2*Pi))*(stor::Ly/stor::Lz)*log10(1 + (stor::Lz/stor::Ly)*(stor::Lz/stor::Ly));
        stor::Nz = 1 - (2/Pi)*atan(stor::Lz/stor::Ly) + (1/(2*Pi))*(stor::Ly/stor::Lz)*log10(1 + (stor::Lz/stor::Ly)*(stor::Lz/stor::Ly)) - (1/(2*Pi))*(stor::Lz/stor::Ly)*log10(1 + (stor::Ly/stor::Lz)*(stor::Ly/stor::Lz));
    }

    // Recompute every quantity that depends on Ms (including the integration
    // prefactors). Must be called whenever Ms changes.
    void update_magnetic_parameters(){
        stor::muMs = stor::mu0*stor::Ms;
        stor::H_demag = stor::Ms*(stor::Nz - stor::Ny);
        calculate::prefac1 =(-stor::alpha*stor::gamma)/((1+pow(stor::alpha,2))*2*stor::Ms*stor::Lz*stor::Ly);
        calculate::prefac2 =stor::mu0*stor::gamma*stor::H_demag/2.0;
        calculate::prefac3 =-stor::gamma/((1+ stor::alpha*stor::alpha)*2*stor::Ms*stor::Lz*stor::Ly);
        calculate::prefac4 =-(stor::gamma*stor::alpha*stor::mu0*stor::H_demag)/(2+2*stor::alpha*stor::alpha);
    }

    // Temperature at time t (in s)
    double T_of_t(double time){
        if( !stor::use_T_pulse ) return stor::T0;
        const double dt = time - stor::T_pulse_t0;
        return stor::T0 + stor::T_pulse_peak*exp(-dt*dt/(2.0*stor::T_pulse_sigma*stor::T_pulse_sigma));
    }

    // Set the current temperature and update Ms and everything depending on it
    void set_temperature(double T){
        stor::T_sim = T;
        stor::T_e = T;
        stor::T_l = T;
        if( stor::use_Ms_T ){
            stor::Ms = stor::Ms_T(T);
            stor::update_magnetic_parameters();
        }
    }

    // Called by the integrators at every evaluation of the right-hand side,
    // i.e. also at the intermediate Runge-Kutta stages.
    void update_temperature(double time){
        if( stor::use_T_pulse ){
            stor::set_temperature(stor::T_of_t(time));
        }
        else if( ttm::enabled ){
            ttm::evaluate(time);
            stor::T_e = ttm::Te;
            stor::T_l = ttm::Tl;
            // Spin temperature: drives Ms(T) and the thermal noise
            stor::T_sim = ttm::spin_follows_electrons ? ttm::Te : ttm::Tl;
            if( stor::use_Ms_T ){
                stor::Ms = ( ttm::tau_m > 0.0 ) ? ttm::m*stor::Ms0 : stor::Ms_T(stor::T_sim);
                if( stor::Ms < 1e-3*stor::Ms0 ){
                    std::cerr << "Error: at t = " << time << " s the spin temperature is " << stor::T_sim
                              << " K and Ms = " << stor::Ms << " A/m (< 0.1% of Ms0, Tc = " << stor::Tc << " K)."
                              << " The 1D DW model is not defined for a demagnetised wire."
                              << " Reduce the pulse energy or use a finite TTM_tau_m." << std::endl;
                    exit(-1);
                }
                stor::update_magnetic_parameters();
            }
        }
    }

    // Append the temperatures and Ms at the given time to an output line,
    // keeping the stream formatting unchanged for the other columns.
    void write_temperature_columns(std::ostream &out, double time){
        stor::update_temperature(time);
        std::ios::fmtflags flags = out.flags();
        std::streamsize prec = out.precision();
        out << std::scientific << std::setprecision(8)
            << "\t" << stor::T_e << "\t" << stor::T_l << "\t" << stor::Ms;
        out.flags(flags);
        out.precision(prec);
    }

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

        // Temperature dependent saturation magnetisation
        // Ms from the input file is interpreted as Ms0 (T = 0 K) when use_Ms_T=1
        stor::Ms0 = stor::Ms;
        stor::T0 = stor::T_sim;
        stor::use_Ms_T = ( inputs.get<int>("use_Ms_T", 0) > 0 );
        if( stor::use_Ms_T ){
            stor::Tc = inputs.get<double>("Tc", 0.0);
            stor::Ms_exponent = inputs.get<double>("Ms_exponent", 0.0);
            if( stor::Tc <= 0.0 || stor::Ms_exponent <= 0.0 ){
                std::cerr << "Error: use_Ms_T=1 requires Tc>0 and Ms_exponent>0 in the input file." << std::endl;
                exit(-1);
            }
        }

        // Gaussian temperature pulse
        stor::use_T_pulse = ( inputs.get<int>("use_T_pulse", 0) > 0 );
        if( stor::use_T_pulse ){
            stor::T_pulse_peak  = inputs.get<double>("T_pulse_peak", 0.0);
            stor::T_pulse_t0    = inputs.get<double>("T_pulse_t0", 0.0);
            stor::T_pulse_sigma = inputs.get<double>("T_pulse_sigma", 0.0);
            if( stor::T_pulse_sigma <= 0.0 ){
                std::cerr << "Error: use_T_pulse=1 requires T_pulse_sigma>0." << std::endl;
                exit(-1);
            }
        }

        // The maximum temperature reached must stay below Tc
        const double T_max = stor::T0 + ( stor::use_T_pulse ? std::max(0.0, stor::T_pulse_peak) : 0.0 );
        if( stor::use_Ms_T && T_max >= stor::Tc ){
            std::cerr << "Error: maximum temperature (" << T_max << " K) >= Tc (" << stor::Tc
                      << " K); Ms vanishes and the 1D DW model is not defined." << std::endl;
            exit(-1);
        }

        stor::compute_demag_factors();
        stor::update_magnetic_parameters();
        stor::freq = inputs.get<double>("f");
        stor::P= inputs.get<double>("P");
        stor::j_dens= inputs.get<double>("j");
        stor::beta= inputs.get<double>("beta");
        stor::omega=2*Pi*stor::freq;
        integrate::Dt = inputs.get<double>("Dt");
        integrate::totaltime = inputs.get<double>("totaltime");
        integrate::out_time = inputs.get<double>("out_time");
        integrate::scheme = inputs.get<std::string>("integrator", "RK4");
        stor::program = inputs.get<std::string>("program", "Benchmark2");

        // Two-temperature model for ultrafast laser heating
        ttm::enabled = ( inputs.get<int>("use_TTM", 0) > 0 );
        if( ttm::enabled ){
            if( stor::use_T_pulse ){
                std::cerr << "Error: use_TTM=1 and use_T_pulse=1 are mutually exclusive." << std::endl;
                exit(-1);
            }
            ttm::gamma_e  = inputs.get<double>("TTM_gamma_e", 0.0);
            ttm::C_l      = inputs.get<double>("TTM_C_l", 0.0);
            ttm::G        = inputs.get<double>("TTM_G", 0.0);
            ttm::t0       = inputs.get<double>("TTM_t0", 0.0);
            ttm::sigma    = inputs.get<double>("TTM_sigma", 0.0);
            ttm::tau_cool = inputs.get<double>("TTM_tau_cool", 0.0);
            ttm::tau_m    = inputs.get<double>("TTM_tau_m", 0.0);
            ttm::dt_max   = inputs.get<double>("TTM_dt_max", 1e-14);
            std::string spin_T = inputs.get<std::string>("TTM_spin_temperature", "electron");
            if( spin_T == "electron" ) ttm::spin_follows_electrons = true;
            else if( spin_T == "phonon" ) ttm::spin_follows_electrons = false;
            else {
                std::cerr << "Error: TTM_spin_temperature must be 'electron' or 'phonon'." << std::endl;
                exit(-1);
            }
            // Pulse energy: either the peak absorbed power density P0 (W/m^3) or the
            // absorbed fluence F (J/m^2), assumed uniformly absorbed over the thickness Lz:
            //   P0 = F / (Lz * sigma * sqrt(2*pi))
            const double fluence = inputs.get<double>("TTM_fluence", 0.0);
            ttm::P0 = inputs.get<double>("TTM_P0", 0.0);
            if( fluence > 0.0 && ttm::P0 > 0.0 ){
                std::cerr << "Error: give either TTM_P0 or TTM_fluence, not both." << std::endl;
                exit(-1);
            }
            if( ttm::sigma <= 0.0 || ttm::gamma_e <= 0.0 || ttm::C_l <= 0.0 || ttm::G <= 0.0 ){
                std::cerr << "Error: use_TTM=1 requires TTM_gamma_e, TTM_C_l, TTM_G and TTM_sigma > 0." << std::endl;
                exit(-1);
            }
            if( stor::T0 <= 0.0 ){
                std::cerr << "Error: use_TTM=1 requires Temperature > 0 (electron heat capacity gamma_e*Te)." << std::endl;
                exit(-1);
            }
            if( ttm::tau_m > 0.0 && !stor::use_Ms_T ){
                std::cerr << "Error: TTM_tau_m > 0 requires use_Ms_T=1." << std::endl;
                exit(-1);
            }
            if( fluence > 0.0 ) ttm::P0 = fluence/(stor::Lz*ttm::sigma*sqrt(2.0*Pi));
            ttm::initialise( stor::T0 );
        }
        stor::set_temperature(stor::T0);
        stor::update_temperature(0.0);

        stor::Nwires = inputs.get<int>("Nwires");
	stor::x_coord.assign(stor::Nwires, 8.70547e-08);
        stor::phi_coord.assign(stor::Nwires, 0.0);
        stor::V0_mdw.assign(stor::Nwires, stor::V0);
        stor::H_DW.assign(stor::Nwires, 0.0);
        stor::H_dd.assign(stor::Nwires, 0.0);

        stor::H_const = inputs.get<double>("H_const", 0.0);

        stor::use_DW_coupling = false;
        if ( inputs.get<int>("use_DW_coupling", 0) > 0)
            stor::use_DW_coupling = true;
        std::cout << "Use DW coupling = " << stor::use_DW_coupling << std::endl;


        stor::x_coord[0] = 10e-9;
        stor::x_coord[1] = 10.0e-9;
        //stor::x_coord[2] = 0.0;

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



        if( inputs.get<int>("use_edge_roughness", 0) > 0){
            stor::use_edge = true;
            stor::H_edge_max = inputs.get<double>("H_edge_max", 0.0);
            stor::edge_scale = inputs.get<double>("edge_scale", 20e-9); // discretisation of the edge roughness
            int Nx = std::round( stor::L / stor::edge_scale);
            std::cerr << stor::edge_scale << "  " << Nx << std::endl;
            stor::H_edge.assign(stor::Nwires, Nx, 0.0);
            stor::E_edge.assign(stor::Nwires, Nx+1, 0.0);

            int edge_seed = inputs.get<int>("edge_seed", 0);

            // if a seed if set for the edge roughness we reset the RNG
            if( edge_seed > 0)
                calculate::seed_rng(edge_seed);


            for( int i = 0; i < stor::Nwires; i++){
                stor::E_edge(i,0) = calculate::Uniform();
                    std::cerr << i << "  " << -1 << "  " << stor::E_edge(i,0) << "  " << 0.0 << std::endl;
                for( int j = 0; j < Nx; j++){
                    stor::E_edge(i, j+1) = calculate::Uniform();
                    stor::H_edge(i,j) = (E_edge(i,j) - E_edge(i,j+1))*stor::H_edge_max;
                    std::cerr << i << "  " << j << "  " << stor::E_edge(i,j+1) << "  " << stor::H_edge(i,j) << std::endl;
                }
            }

            // After the edge roughness has been calculated reset the RNG to
            // its originial random seed
            if( edge_seed > 0)
                calculate::seed_rng(calculate::random_seed);
        }



        if( stor::Nwires > 1 ){
        integrate::multi_dw::setup(stor::Nwires);
	    stor::V0_mdw[1]=stor::V0;
	    }

        std::cout<<"The program has been initialized with following parameters:"<<std::endl;
        if( stor::use_Ms_T ){
            std::cout<<"Ms(T) = Ms0*(1-T/Tc)^"<<stor::Ms_exponent<<", Ms0 = "<<stor::Ms0<<" A/m, Tc = "<<stor::Tc<<" K"<<std::endl;
            std::cout<<"Saturation, Ms(T = "<<stor::T_sim<<" K) = "<<stor::Ms<<" A/m (at t = 0)"<<std::endl;
        }
        if( stor::use_T_pulse ){
            std::cout<<"Gaussian T pulse: T0 = "<<stor::T0<<" K, peak rise = "<<stor::T_pulse_peak
                     <<" K, t0 = "<<stor::T_pulse_t0<<" s, sigma = "<<stor::T_pulse_sigma<<" s"<<std::endl;
            if( stor::use_Ms_T )
                std::cout<<"Ms at pulse maximum = "<<stor::Ms_T(stor::T0+stor::T_pulse_peak)<<" A/m"<<std::endl;
        }
        if( ttm::enabled ){
            std::cout<<"Two-temperature model: gamma_e = "<<ttm::gamma_e<<" J/m^3/K^2, C_l = "<<ttm::C_l
                     <<" J/m^3/K, G = "<<ttm::G<<" W/m^3/K"<<std::endl;
            std::cout<<"  Laser: P0 = "<<ttm::P0<<" W/m^3, t0 = "<<ttm::t0<<" s, sigma = "<<ttm::sigma
                     <<" s, absorbed energy density = "<<ttm::P0*ttm::sigma*sqrt(2.0*Pi)<<" J/m^3"<<std::endl;
            std::cout<<"  tau_cool = "<<ttm::tau_cool<<" s, tau_m = "<<ttm::tau_m<<" s, spin temperature = "
                     <<( ttm::spin_follows_electrons ? "electron" : "phonon" )<<std::endl;
        }
        else
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
        std::cout<<"Predicted critical current = " << stor::gamma*stor::mu0*stor::H_demag*0.5*calculate::calculate_DW(0.78539816339744830962) << std::endl;
        std::cout<<"=====================================================================<"<<std::endl;

        // Initialise the parameters for the integration
        calculate::zeeman_prefac1 = stor::gamma*stor::mu0*stor::alpha/(stor::alpha*stor::alpha+1.0);
        calculate::zeeman_prefac2 = stor::gamma*stor::mu0/(1.0 + stor::alpha*stor::alpha);
    }
}//end of namespace
