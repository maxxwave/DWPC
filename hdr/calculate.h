#ifndef __CALCULATE_H__
#define __CALCULATE_H__

#include <vector>

namespace calculate{
    extern unsigned random_seed;
    void seed_rng(const unsigned);
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
	extern double noise( double T, double DW );
	extern double current_prefac;
	extern double current(double time);
	extern double prefac1, prefac2, prefac3, prefac4, zeeman_prefac1, zeeman_prefac2;
	double x_t(double DWs, double phi, double phi_t);
	double phi_t(double dEx, double phi_rk, double H);
	extern double DW_coupling(std::vector <double> &);
    void gradient( double &, double &, double, double, const double);
    void gradient( std::vector<double> &, std::vector<double> &, std::vector<double> &, std::vector<double> &, const double);
    void BDF( double &, double);
    void BDF1( double &, double);
    void BDF2( double &, double);
    void BDF3( double &, double);
    void noise_gradient( double&, double&, double, double, double, double);
    double Normal();
}
#endif
