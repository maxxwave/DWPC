// Author: Razvan Ababei	
// University of Sheffield
// date 23rd of August, 2019
//
//
// This file is dedicated to implement the basic Euler integration scheme
// for DW propagation code. This file comes without any warranty
//
//

#include <iostream>
#include <vector>
#include <cmath>
#include "../hdr/storage.hdr"
#include "../hdr/create.hdr"
#include "../hdr/euler_integrator.hdr"
//#include "../hdr/update.hdr"


namespace integrate{
	
	const double one_rad=0.0175;	
	// defining local variables
	double x_euler=0.0;
	double phi_euler=0.0;	
	double vx_euler=0.0;
	double phi_dt_euler=0.0;

	// setting the time step of integration
	double Dt=1e-15; // in seconds

	// Defining some prefactors where we incorporate the constants in order to not be called each time in the loop
	const double prefac1=(-stor::alpha*stor::gamma)/
                        ((1+pow(stor::alpha,2))*2*stor::Ms*stor::Lz*stor::Ly);
	const double prefac2=stor::mu0*stor::gamma*stor::H_demag/(2.0+2.0*stor::alpha*stor::alpha);
	const double prefac3=(-stor::gamma)/
                        ((1+pow(stor::alpha,2))*2*stor::Ms*stor::Lz*stor::Ly);
	const double prefac4=-stor::gamma*stor::alpha*stor::mu0*stor::H_demag/(2*(1+stor::alpha*stor::alpha));
	const double zeeman_prefac1= stor::gamma*stor::mu0*stor::alpha/(stor::alpha*stor::alpha+1.0);
	const double zeeman_prefac2=stor::gamma*stor::mu0/(1.0 + stor::alpha*stor::alpha);	
		
	double euler(){
		
		//std::cout<<prefac1<<"\t"<<prefac2<<"\t"<<prefac3<<"\t"<<prefac4<<std::endl;

		// transfer the coordinates into eurler variables
		x_euler=stor::x_dw;
		phi_euler=stor::phi_dw;
		//std::cout<<x_euler<<"\t"<<stor::x_dw<<std::endl;
		
		// calculate the domain width 
		create::update_energy(x_euler);
		create::calculate_DW(phi_euler);	
		//std::cout<<stor::Ex<<"\t"<<stor::dEx<<"\t"<<stor::Dw_size<<std::endl;

		// calculate the speed and angular speed
		vx_euler= prefac1*stor::Dw_size*stor::dEx +
		       	prefac2*stor::Dw_size*sin(2*phi_euler*one_rad) +
		       	zeeman_prefac1*stor::Dw_size*stor::V;

		phi_dt_euler=prefac3*stor::dEx +
		       	prefac4*sin(2*phi_euler*one_rad)+
			zeeman_prefac2*stor::V;
		
		// transfer the speed values into storage
		stor::vx=vx_euler;
		stor::phi_dt=phi_dt_euler;

		//calculate the final step using trivial euler method
		stor::x_dw= x_euler + Dt*vx_euler;
		stor::phi_dw=phi_euler + Dt*phi_dt_euler; 
		
		//std::cout<<stor::x_dw<<"\t"<<stor::Dw_size<<"\t"<<stor::phi_dw<<std::endl;
		//std::cout<<prefac1<<"\t"<<prefac2<<"\t"<<"DW_SIZE"<<stor::Dw_size<<"\t"<<"DEX"<<"\t"<<stor::dEx<<"\t"<<stor::x_dw<<"\t"<<stor::phi_dw<<"\t"<<stor::H_demag<<"\t"<<vx_euler<<"\t"<<phi_dt_euler<<std::endl;

	}// end of euler function


}// end of namespace 