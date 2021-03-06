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

}
#endif
