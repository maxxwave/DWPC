#include <iostream>
#include <vector>
#include "../hdr/storage.h"
#include "../hdr/calculate.h"
#include "../hdr/euler_integrator.h"
#ifndef __RC_H__
#define __RC_H__

namespace reservoir{
	extern std::vector<double> H0;
	extern int no_nodes;
	extern double tau, theta;
	extern double y_p;
	extern const double r;
	extern const double sigma;
	extern std::vector <double> W;
	extern std::vector <double> s_x;
	extern std::vector <double> t_p;
	extern std::vector <double> H_class;
	extern std::vector <double> mask_array;
	extern double sigmoid(double x);
	extern double mask_values();
	double classification(std::vector<double>&, std::vector<double>&);
	double training( std::vector<double>&, std::vector<double>&);
	//extern double V_min, V_max;
	extern long int no_steps_per_node;
	extern double time;
	extern double mask();
	extern double oscillation_response(double Hi);
    extern double field_map();
	void get_input_data(std::string&, std::vector<double> &, std::vector<double> &);
    int run();
}
#endif
