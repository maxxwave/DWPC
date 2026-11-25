// Two-temperature model for ultrafast laser heating in DWPC.
// See hdr/ttm.h for the equations.
//
// The TTM is integrated with its own adaptive RK4 sub-steps, because the
// electron-phonon and pulse time scales (~10-1000 fs) are usually much shorter
// than the DW integration step. The TTM depends only on time, so each request
// is served by integrating forward from the latest stored state at or before
// the requested time. A short history of states covers integrator stages that
// step slightly back in time (e.g. the last RK45 stage); a request earlier than
// the whole history (a program restarting time at 0) re-integrates from t = 0.

#include <algorithm>
#include <deque>
#include <cmath>
#include <iostream>
#include <stdlib.h>
#include "../hdr/storage.h"
#include "../hdr/ttm.h"

namespace ttm{

    bool enabled = false;
    double gamma_e = 0.0, C_l = 0.0, G = 0.0, P0 = 0.0;
    double t0 = 0.0, sigma = 0.0, tau_cool = 0.0, tau_m = 0.0;
    double dt_max = 1e-14;
    bool spin_follows_electrons = true;

    double Te = 0.0, Tl = 0.0, m = 1.0;

    namespace {
        struct state_t { double t, Te, Tl, m; };
        state_t initial;
        std::deque<state_t> history;       // increasing in time
        const std::size_t history_size = 32;
        double T_bath = 0.0;

        void rhs( const double t, const double y[3], double dy[3]){
            const double Ce = gamma_e*y[0];
            dy[0] = ( -G*(y[0]-y[1]) + power(t) ) / Ce;
            dy[1] = G*(y[0]-y[1]) / C_l;
            if( tau_cool > 0.0 ) dy[1] -= (y[1] - T_bath)/tau_cool;
            if( tau_m > 0.0 ){
                const double Ts = spin_follows_electrons ? y[0] : y[1];
                const double m_eq = stor::Ms_T(Ts)/stor::Ms0;
                dy[2] = (m_eq - y[2])/tau_m;
            } else dy[2] = 0.0;
        }

        // Step size limited by the fastest relaxation time and the pulse width
        double step_size( const state_t &s){
            double h = dt_max;
            h = std::min(h, 0.02*gamma_e*s.Te/G);
            h = std::min(h, 0.02*C_l/G);
            if( tau_cool > 0.0 ) h = std::min(h, 0.02*tau_cool);
            if( tau_m > 0.0 )    h = std::min(h, 0.02*tau_m);
            if( std::fabs(s.t - t0) < 8.0*sigma ) h = std::min(h, sigma/20.0);
            return h;
        }

        void advance( state_t &s, const double t_target){
            while( s.t < t_target ){
                const double h = std::min( step_size(s), t_target - s.t );
                double y[3] = { s.Te, s.Tl, s.m }, yt[3];
                double k1[3], k2[3], k3[3], k4[3];
                rhs(s.t, y, k1);
                for(int i=0;i<3;i++) yt[i] = y[i] + 0.5*h*k1[i];
                rhs(s.t + 0.5*h, yt, k2);
                for(int i=0;i<3;i++) yt[i] = y[i] + 0.5*h*k2[i];
                rhs(s.t + 0.5*h, yt, k3);
                for(int i=0;i<3;i++) yt[i] = y[i] + h*k3[i];
                rhs(s.t + h, yt, k4);
                s.Te += h*(k1[0] + 2*k2[0] + 2*k3[0] + k4[0])/6.0;
                s.Tl += h*(k1[1] + 2*k2[1] + 2*k3[1] + k4[1])/6.0;
                s.m  += h*(k1[2] + 2*k2[2] + 2*k3[2] + k4[2])/6.0;
                s.t  += h;
                if( !(s.Te > 0.0) || !std::isfinite(s.Te) || !std::isfinite(s.Tl) ){
                    std::cerr << "Error: TTM integration became unstable at t = " << s.t
                              << " s (Te = " << s.Te << ", Tl = " << s.Tl << "). Reduce TTM_dt_max." << std::endl;
                    exit(-1);
                }
            }
        }
    }

    double power( double time){
        const double dt = time - t0;
        return P0*exp(-dt*dt/(2.0*sigma*sigma));
    }

    void initialise( double T0_bath){
        T_bath = T0_bath;
        double m0 = 1.0;
        if( stor::use_Ms_T ) m0 = stor::Ms_T(T0_bath)/stor::Ms0;
        initial.t = 0.0; initial.Te = T0_bath; initial.Tl = T0_bath; initial.m = m0;
        history.clear();
        history.push_back(initial);
        Te = T0_bath; Tl = T0_bath; m = m0;
    }

    void evaluate( double time){
        state_t s;
        if( time >= history.back().t ){
            s = history.back();
            advance( s, time );
            if( s.t > history.back().t ){
                history.push_back(s);
                if( history.size() > history_size ) history.pop_front();
            }
        }
        else {
            // latest stored state at or before the requested time
            int i = int(history.size()) - 1;
            while( i >= 0 && history[i].t > time ) i--;
            if( i < 0 ){
                history.clear();
                history.push_back(initial);
                s = initial;
                advance( s, time );
                history.push_back(s);
            } else {
                s = history[i];
                advance( s, time );
            }
        }
        Te = s.Te; Tl = s.Tl; m = s.m;
    }
}
