// Author: Razvan Ababei	
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
#include "../hdr/create.hdr"
#include "../hdr/euler_integrator.hdr"
#include "../hdr/runge_kutta.hdr"

namespace integrate{


	// in this function we implemented the 4th order range kutta integration scheme
	// yn+1=yn + 1/6(k1 +2k2 + 2k3 +k4)
	// k1 = h f(tn, yn)
	// k2 = h f(tn + h/2, yn + k1/2)
	// k3 = h f(tn + h/2, yn + k2/2)
	// k4 = h f(tn + h, yn + k3)
	// h = xn+1-xn
	//
	// we define two function for speed and angular speed
	double phi_t(double dEx, double phi_rk, double H){
		return prefac3*dEx+prefac4*sin(2*phi_rk)+zeeman_prefac2*V;	
	}
	double x_t(double DWs, double phi){
		return 
	}
	double runge_kutta(){
	}


}// end of namespace
