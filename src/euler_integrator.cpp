// Author: Razvan Ababei
// University of Sheffield
// date 23rd of August, 2019
//
//
// This file is dedicated to implement the basic Euler integration scheme
// for DW propagation code. This file comes without any warranty
//
// This version is in agreement with Boulle et al.
//

#include <iostream>
#include <vector>
#include <cmath>
#include "../hdr/storage.h"
#include "../hdr/calculate.h"
#include "../hdr/euler_integrator.h"
//#include "../hdr/update.hdr"

using namespace calculate;
namespace integrate{

	// defining local variables
	double x_euler=0.0;
	double phi_euler=0.0;
	double vx_euler=0.0;
	double phi_dt_euler=0.0;

	// setting the time step of integration
	double Dt=1e-15; // in seconds
	double totaltime = 1e-9;
	double out_time = 1e-10;
	std::string scheme;


	double euler( double &time){

		//calculate::Zeeman(time);
		//std::cout<<prefac1<<"\t"<<prefac2<<"\t"<<prefac3<<"\t"<<prefac4<<std::endl;

		// transfer the coordinates into eurler variables
		x_euler=stor::x_dw;
		phi_euler=stor::phi_dw;
		//std::cout<<x_euler<<"\t"<<stor::x_dw<<std::endl;

		// calculate the domain width
		calculate::update_energy_antinotches(x_euler);
		calculate::calculate_DW(phi_euler);

		//std::cout<<stor::Ex<<"\t"<<stor::dEx<<"\t"<<stor::Dw_size<<std::endl;

		// calculate the speed and angular speed
		phi_dt_euler = prefac3*stor::dEx +
		       	       prefac4*sin(2*phi_euler)+
			       zeeman_prefac2*stor::V;

		//vx_euler = stor::Dw_size*stor::gamma*stor::mu0*stor::V/stor::alpha+
			   //stor::Dw_size*phi_dt_euler/stor::alpha;

		vx_euler = prefac2*sin(2*phi_euler)*stor::Dw_size + stor::alpha*stor::Dw_size*phi_dt_euler;
		//vx_euler= prefac1*stor::dEx + prefac2*stor::Dw_size+ zeeman_prefac1*stor::V;

		// transfer the speed values into storage
		stor::vx=vx_euler;
		stor::phi_dt=phi_dt_euler;

		//calculate the final step using trivial euler method
		stor::x_dw= x_euler + Dt*vx_euler;
		stor::phi_dw=phi_euler + Dt*phi_dt_euler;
		//std::cout<<stor::x_dw<<"\t"<<stor::Dw_size<<"\t"<<stor::phi_dw<<std::endl;
		/*std::cout<<prefac1<<"\t"<<prefac2<<"\t"<<
			"DW_SIZE"<<"\t"<<stor::Dw_size<<"\t"
			<<"DEX"<<"\t"<<stor::dEx<<"\t"
			<<stor::x_dw<<"\t"<<stor::phi_dw<<"\t"
			<<stor::H_demag<<"\t"<<vx_euler<<"\t"<<phi_dt_euler<<"\t"
			<<prefac1*stor::Dw_size*stor::dEx<<"\t"<<
			prefac2*stor::Dw_size*sin(2*phi_euler*one_rad)<<"\t"<<
			zeeman_prefac1*stor::Dw_size*stor::V<<"\t"<<
			std::endl;
		*/
		time += integrate::Dt;
	return 1;
	}// end of euler function


}// end of namespace
