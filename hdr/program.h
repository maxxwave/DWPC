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
}
#endif

