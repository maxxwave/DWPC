// Utility file
// Created:     2/2/2023
// Modified:    2/2/2023
// Author:      Matthew Ellis
//

#ifndef __UTILITY_H__
#define __UTILITY_H__

#include <cassert>

template < typename T >
inline T sqr ( const T &x ) { return x * x; }

template < typename T >
inline T cube ( const T &x) { return x * x * x; }

class double3 {
    public:
        double _data[3];

        double3 () { _data[0]=0.0; _data[1]=0.0; _data[2]=0.0;}
        double3 (double val) { _data[0]=val; _data[1]=val; _data[2]=val;}
        double3 (double a, double b, double c) { _data[0]=a; _data[1]=b; _data[2]=b;}
        ~double3 () {};

        double&  operator[] ( int i ) {
            assert( (i < 3) && ( i >= 0) );
        return _data[i];
        }
};

namespace util {

    const double Pi =  3.1415926535897932384626;
    const double Pi_sqr =  Pi*Pi;
    const double two_Pi = 2*Pi;
    const double root_Pi = sqrt(Pi);

    // Gyromagnetic ratio in per ns per T
    const double gyro_ns = 176.0;

    const double mu_0 = 1.256637061e-6;

    const double kB = 1.38e-23;

    const double e = 1.60217653e-19;

} //end of namespace util
#endif

