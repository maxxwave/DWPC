// DWPC - Domain Wall Propagation Code
// Copyright (C) 2019-2026 Razvan-Vasile Ababei and Matt Ellis
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This program is free software: you can redistribute it and/or modify it under the terms of
// the GNU General Public License as published by the Free Software Foundation, either version 3
// of the License, or (at your option) any later version. It is distributed WITHOUT ANY WARRANTY;
// see the file LICENSE for details. If you use DWPC in published work, please cite it as
// described in CITATION.cff.

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
	extern void mask_values();
	double classification(std::vector<double>&, std::vector<double>&);
	double training( std::vector<double>&, std::vector<double>&);
	//extern double V_min, V_max;
	extern long int no_steps_per_node;
	extern double time;
	extern void oscillation_response(double Hi);
	extern double field_map();
	void get_input_data(std::string&, std::vector<double> &, std::vector<double> &);
	int run();
	int run_spoken_recognition();
    int run_field_sequence();
    int run_transient();
    int multi_wires_dynamics();
  
}
#endif
