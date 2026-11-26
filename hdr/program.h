// DWPC - Domain Wall Propagation Code
// Copyright (C) 2019-2026 Razvan-Vasile Ababei and Matt Ellis
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This program is free software: you can redistribute it and/or modify it under the terms of
// the GNU General Public License as published by the Free Software Foundation, either version 3
// of the License, or (at your option) any later version. It is distributed WITHOUT ANY WARRANTY;
// see the file LICENSE for details. If you use DWPC in published work, please cite it as
// described in CITATION.cff.

#ifndef __PROGRAMS_H__
#define __PROGRAMS_H__

#include <iostream>
#include "../hdr/calculate.h"
#include "../hdr/storage.h"
#include "../hdr/euler_integrator.h"

namespace programs{
	extern double bifurcation();
	extern double benchmark1();
	extern double benchmark2();
	extern double benchmark3();
	extern double show_potential();
	extern double spin_current1();
	extern double Lyapunov_exponent(); 
}
#endif

