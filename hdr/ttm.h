#ifndef __TTM_H__
#define __TTM_H__

#include <string>

// Two-temperature model (TTM) for ultrafast laser heating.
//
//   gamma_e*Te * dTe/dt = -G*(Te - Tl) + P(t)
//   C_l        * dTl/dt =  G*(Te - Tl) - C_l*(Tl - T0)/tau_cool
//   P(t) = P0*exp(-(t - t0)^2 / (2*sigma^2))
//
// Optional finite magnetisation response (tau_m > 0):
//   dm/dt = (m_eq(Ts) - m)/tau_m,  m_eq(T) = Ms(T)/Ms0,  Ms = m*Ms0
// where Ts is the spin temperature (electron or phonon temperature).
namespace ttm{

    extern bool enabled;
    extern double gamma_e;   // J m^-3 K^-2, electron heat capacity Ce = gamma_e*Te
    extern double C_l;       // J m^-3 K^-1, lattice heat capacity
    extern double G;         // W m^-3 K^-1, electron-phonon coupling
    extern double P0;        // W m^-3, peak absorbed power density
    extern double t0;        // s, time of the pulse maximum
    extern double sigma;     // s, standard deviation of the Gaussian pulse
    extern double tau_cool;  // s, lattice cooling to the substrate (0 = off)
    extern double tau_m;     // s, magnetisation relaxation time (0 = instantaneous)
    extern double dt_max;    // s, maximum internal TTM step
    extern bool spin_follows_electrons;

    // State evaluated at the last requested time
    extern double Te, Tl, m;

    void initialise(double T0_bath);
    void evaluate(double time);   // sets Te, Tl, m at the given time
    double power(double time);
}
#endif
