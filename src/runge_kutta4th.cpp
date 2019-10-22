// Author: Razvan Ababei & Matt Ellis	
// University of Sheffield
// date 21st of August, 2019
//
// This file is dedicated to implement the integrator based on Runge Kutta method 
//

#include <iostream>
#include <vector>
#include <cmath>

#include "../hdr/update.hdr"
#include "../hdr/storage.hdr"
#include "../hdr/euler_integrator.hdr"
#include "../hdr/calculate.hdr"
namespace integrate{


	// in this function we implemented the 4th order range kutta integration scheme
	// yn+1=yn + 1/6(k1 +2k2 + 2k3 +k4)
	// k1 = h f(tn, yn)
	// k2 = h f(tn + h/2, yn + k1/2)
	// k3 = h f(tn + h/2, yn + k2/2)
	// k4 = h f(tn + h, yn + k3)
	// h = xn+1-xn
	
	double k1=0.0, k2=0.0, k3=0.0, k4=0.0;
	double kp1=0.0, kp2=0.0, kp3=0.0, kp4=0.0;
	double phi_k=0.0, x_k=0.0;


	// we define two function for speed and angular speed
	double phi_t(double dEx, double phi_rk, double H){
		return prefac3*dEx+prefac4*sin(2*phi_rk)+ zeeman_prefac2*H;	
	}
	double x_t(double DWs, double phi, double phi_t){
		return prefac2*sin(2*phi)*DWs + stor::alpha*DWs*phi_t; 
	}
	double runge_kutta(){
		phi_k=stor::phi_dw;
		x_k=stor::x_dw;

		// calculate the domain width 
                calculate::update_energy_antinotches(x_k);
                calculate::calculate_DW(phi_k);

		kp1=integrate::Dt*phi_t(stor::dEx, phi_k, stor::V);
		k1=integrate::Dt*x_t(stor::Dw_size, phi_k, phi_t(stor::dEx, phi_k, stor::V));
		calculate::update_energy_antinotches(x_k+0.5*k1);
                calculate::calculate_DW(phi_k+0.5*kp1);

		kp2=integrate::Dt*phi_t(stor::dEx, phi_k+0.5*kp1, stor::V);
		k2=integrate::Dt*x_t(stor::Dw_size,phi_k+0.5*kp1, phi_t(stor::dEx, phi_k, stor::V));
		calculate::update_energy_antinotches(x_k+0.5*k2);
                calculate::calculate_DW(phi_k+0.5*kp2);

		kp3=integrate::Dt*phi_t(stor::dEx, phi_k+0.5*kp2, stor::V);
		k3=integrate::Dt*x_t(stor::Dw_size,phi_k+0.5*kp2,phi_t(stor::dEx, phi_k+0.5*kp2, stor::V));
		calculate::update_energy_antinotches(x_k+0.5*k3);
                calculate::calculate_DW(phi_k+0.5*kp3);


		kp4=integrate::Dt*phi_t(stor::dEx, phi_k+kp3, stor::V);
		k4=integrate::Dt*x_t(stor::Dw_size,phi_k+kp3,phi_t(stor::dEx, phi_k+kp3, stor::V));
		stor::vx=k4;
		stor::phi_dt=kp4;
		stor::phi_dw = phi_k + (1.0/6.0)*(kp1 +2*kp2+2*kp3+kp4);
		stor::x_dw = x_k + (1.0/6.0)*(k1 +2*k2+2*k3+k4);

	}


}// end of namespace
