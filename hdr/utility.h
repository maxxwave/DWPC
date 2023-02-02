// Utility file
// Created:     2/2/2023
// Modified:    2/2/2023
// Author:      Matthew Ellis
//

#ifndef __UTILITY_H__
#define __UTILITY_H__

template < typename T >
inline T sqr ( const T &x ) { return x * x; }

template < typename T >
inline T cube ( const T &x) { return x * x * x; }

namespace util {

    const double Pi =  3.1415926535897932384626;
    const double Pi_sqr =  Pi*Pi;
    const double two_Pi = 2*Pi;
    const double root_Pi = sqrt(pi);

    // Gyromagnetic ratio in per ns per T
    const double gyro_ns = 176.0;

    const double mu_0 = 1.256637061e-6;

    const double kB = 1.38e-23;

    const double e = 1.60217653e-19;

} //end of namespace util
#endif

