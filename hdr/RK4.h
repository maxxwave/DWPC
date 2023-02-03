// Runge-Kutta 4 integrator
// Created:     3/2/2023
// Modified:    3/2/2023
// Author:      Matthew Ellis
//

#ifndef __RK4_H__
#define __RK4_H__

#include "../hdr/DW_model.h"
#include "../hdr/utility.h"
#include "../hdr/calculate.h"

class RK4_t {
    public:

    double dt;
    double k1, k2, k3, k4;
	double kp1, kp2, kp3, kp4;
	double phi_k, x_k;

    RK4_t () {};
    ~RK4_t () {};

    void init()
    {
        dt = 1e-2;
    }

    void step (double &q, double &phi, double &t, DW_model_t &model);
    void multistep (const int N, double &q, double &phi, double &t, DW_model_t &model);

};

void RK4_t::step (double &q, double &phi, double &t, DW_model_t &model )
{
    const double q0 = q;
    const double phi0 = phi;

    model.gradient( k1, kp1, q, phi, t);

    model.gradient( k2, kp2, q + 0.5*k1*dt, phi + 0.5*kp1*dt, t + 0.5*dt);
    
    model.gradient( k3, kp3, q + 0.5*k2*dt, phi + 0.5*kp2*dt, t + 0.5*dt);

    model.gradient( k4, kp4, q + k3*dt, phi + kp3*dt, t + dt);

    q = q0 + (k1 + 2*k2 + 2*k3 + k4)*dt/6.0;
    phi = phi0 + (kp1 + 2*kp2 + 2*kp3 + kp4)*dt/6.0;
    t += dt;
}

void RK4_t::multistep( const int N, double &q, double &phi, double &t, DW_model_t &model  )
{
    for (int i = 0; i < N; i++)
        step(q, phi, t, model);
}

#endif