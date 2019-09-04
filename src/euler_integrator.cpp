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
                        ((1+pow(stor::alpha,2))*2*stor::muMs*stor::Lz*stor::Ly);
	const double prefac2=stor::gamma*stor::H_demag/2.0;
	const double prefac3=(-stor::gamma)/
                        ((1+pow(stor::alpha,2))*2*stor::muMs*stor::Lz*stor::Ly);
	const double prefac4=-prefac2*stor::alpha;
	
		
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

		//calculate the speed and angular speed
		vx_euler= prefac1*stor::Dw_size*stor::dEx + prefac2*stor::H_demag*sin(2*phi_euler*one_rad);
		phi_dt_euler=prefac3*stor::dEx + prefac4*stor::H_demag*sin(2*phi_euler*one_rad);
		
		//calculate the final step using trivial euler method
		stor::x_dw= x_euler + Dt*vx_euler;
		stor::phi_dw=phi_euler + Dt*phi_dt_euler; 
		
		//std::cout<<stor::x_dw<<"\t"<<stor::Dw_size<<"\t"<<stor::phi_dw<<std::endl;
		std::cout<<prefac1<<"\t"<<prefac2<<"\t"<<"DW_SIZE"<<stor::Dw_size<<"\t"<<"DEX"<<stor::dEx<<"\t"<<stor::x_dw<<"\t"<<stor::phi_dw<<"\t"<<stor::H_demag<<"\t"<<vx_euler<<std::endl;

	}// end of euler function


}// end of namespace 
