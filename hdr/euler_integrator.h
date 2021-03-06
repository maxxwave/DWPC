#include <iostream>
#include <vector>
#include <string>
#include "storage.h"
//#include "create.h"
#ifndef __EULER_H__
#define __EULER_H__

namespace integrate {
	extern std::string scheme;
	extern double totaltime;
	extern double out_time;
	extern double Dt;
	extern double prefac1, prefac2, prefac3, prefac4, zeeman_prefac1, zeeman_prefac2;
	extern double k1, k2, k3, k4;
	extern double kp1, kp2, kp3, kp4;
	extern double phi_k, x_k;
	extern double x_euler;
	extern double phi_euler;
	extern double vx_euler;
	extern double phi_dt_euler;
	extern double euler_init();
	extern double euler(double &time);
	extern double x_t(double DWs, double phi, double phi_t);
	extern double phi_t(double dEx, double phi_rk, double H);
	extern double runge_kutta(double &time);
}
#endif
