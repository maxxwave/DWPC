#include <iostream>
#include "storage.h"

#ifndef __CALCULATE_H__
#define __CALCULATE_H__

namespace calculate{
	extern int N;
	extern int a;
	extern double rad;
	extern double c1,c2;
	extern double create();
	extern double initialize();
	extern double calculate_DW(double phi);
	extern double update_energy(double x);
	extern double update_energy_notches(double x);
	extern double update_energy_antinotches(double x);
	extern double Zeeman(double Dt);

	extern double prefac1, prefac2, prefac3, prefac4, zeeman_prefac1, zeeman_prefac2;
	double x_t(double DWs, double phi, double phi_t);
	double phi_t(double dEx, double phi_rk, double H);
    void gradient( double &, double &, double, double, const double);
}
#endif
