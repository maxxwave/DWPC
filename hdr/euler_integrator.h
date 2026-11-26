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
	extern double k1, k2, k3, k4;
	extern double kp1, kp2, kp3, kp4;
	extern double phi_k, x_k;
	extern double x_euler;
	extern double phi_euler;
	extern double vx_euler;
	extern double phi_dt_euler;
	double euler_init();
	double euler(double &time);
	double runge_kutta(double &time);
    int heun(double &time);

    namespace multi_dw {
        void setup (int);
        double runge_kutta( std::vector<double> &,
                std::vector<double>&,
                double &,
                const double);
    }
    namespace RK45 {
        double runge_kutta_45 (double &time, const double dt);
    }
}
#endif
