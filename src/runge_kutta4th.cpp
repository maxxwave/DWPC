// Author: Razvan Ababei & Matt Ellis
// University of Sheffield
// date 21st of August, 2019
//
// This file is dedicated to implement the integrator based on Runge Kutta method
//

#include <iostream>
#include <vector>
#include <cmath>

#include "../hdr/storage.h"
#include "../hdr/euler_integrator.h"
#include "../hdr/calculate.h"

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


	double runge_kutta(double &time){

		// Store the initial positions
        phi_k = stor::phi_dw;
		x_k   = stor::x_dw;

        calculate::gradient( k1, kp1, x_k, phi_k, time);
        calculate::gradient( k2, kp2, x_k + 0.5*Dt*k1, phi_k + 0.5*Dt*kp1, time + 0.5*Dt);
        calculate::gradient( k3, kp3, x_k + 0.5*Dt*k2, phi_k + 0.5*Dt*kp2, time + 0.5*Dt);
        calculate::gradient( k4, kp4, x_k + Dt*k3, phi_k + Dt*kp3, time + Dt);

        stor::phi_dw = phi_k + (kp1 + 2*(kp2 + kp3) + kp4)*Dt/6.0;
		stor::x_dw = x_k + (k1 + 2*(k2 + k3) + k4)*Dt/6.0;
		time += integrate::Dt;

        calculate::gradient( stor::vx, stor::phi_dt, stor::x_dw, stor::phi_dw, time);

          /*
		kp1 = integrate::Dt * calculate::phi_t(stor::dEx, phi_k, stor::V);
		k1 = integrate::Dt * calculate::x_t(stor::Dw_size, phi_k, calculate::phi_t(stor::dEx, phi_k, stor::V));

		calculate::update_energy_antinotches(x_k+0.5*k1);
        calculate::calculate_DW(phi_k+0.5*kp1);
		calculate::Zeeman(time+integrate::Dt*0.5);

		kp2=integrate::Dt*calculate::phi_t(stor::dEx, phi_k+0.5*kp1, stor::V);
		k2=integrate::Dt*calculate::x_t(stor::Dw_size, phi_k+0.5*kp1, calculate::phi_t(stor::dEx, phi_k+0.5*kp1, stor::V));

		calculate::update_energy_antinotches(x_k+0.5*k2);
        calculate::calculate_DW(phi_k+0.5*kp2);
		calculate::Zeeman(time+integrate::Dt*0.5);

		kp3=integrate::Dt*calculate::phi_t(stor::dEx, phi_k+0.5*kp2, stor::V);
		k3=integrate::Dt*calculate::x_t(stor::Dw_size,phi_k+0.5*kp2,calculate::phi_t(stor::dEx, phi_k+0.5*kp2, stor::V));

		calculate::update_energy_antinotches(x_k+k3);
        calculate::calculate_DW(phi_k+kp3);
		calculate::Zeeman(time+integrate::Dt);

		kp4=integrate::Dt*calculate::phi_t(stor::dEx, phi_k+kp3, stor::V);
		k4=integrate::Dt*calculate::x_t(stor::Dw_size,phi_k+kp3,calculate::phi_t(stor::dEx, phi_k+kp3, stor::V));

		stor::vx=calculate::x_t(stor::Dw_size, phi_k, calculate::phi_t(stor::dEx, phi_k, stor::V));
		stor::phi_dt=kp4;
		stor::phi_dw = phi_k + (1.0/6.0)*(kp1 +2*kp2+2*kp3+kp4);
		stor::x_dw = x_k + (1.0/6.0)*(k1 +2*k2+2*k3+k4);
		time += integrate::Dt;
        */

	return 1;
	}// end of function runge_kutta


}// end of namespace
