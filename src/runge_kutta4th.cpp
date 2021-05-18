// Author: Razvan Ababei & Matt Ellis
// University of Sheffield
// date 21st of August, 2019
//
// This file is dedicated to implement the integrator based on Runge Kutta method
//

#include <iostream>
#include <vector>
#include <cmath>

#include "../hdr/storage.h"
#include "../hdr/euler_integrator.h"
#include "../hdr/calculate.h"

namespace integrate{


	// in this function we implemented the 4th order range kutta integration scheme
	// yn+1=yn + 1/6(k1 +2k2 + 2k3 +k4)
	// k1 = h f(tn, yn)
	// k2 = h f(tn + h/2, yn + k1/2)
	// k3 = h f(tn + h/2, yn + k2/2)
	// k4 = h f(tn + h, yn + k3)
	// h = xn+1-xn

	double k1=0.0, k2=0.0, k3=0.0, k4=0.0;
	double kp1=0.0, kp2=0.0, kp3=0.0, kp4=0.0;
	double phi_k=0.0, x_k=0.0;

    double n_x = 0.0, n_phi = 0.0;
    double g1 =0.0, gp1 = 0.0;
    double g2 =0.0, gp2 = 0.0;
    double g3 =0.0, gp3 = 0.0;
    double g4 =0.0, gp4 = 0.0;


    int heun(double &time){


		// Store the initial positions
        phi_k = stor::phi_dw;
        x_k   = stor::x_dw;


        n_x = calculate::Normal();
        n_phi = calculate::Normal();


        calculate::gradient( k1, kp1, x_k, phi_k, time);
        calculate::noise_gradient( g1, gp1, x_k, phi_k, n_x, n_phi);

        calculate::gradient( k2, kp2, x_k + Dt*(k1+g1), phi_k + Dt*(kp1+gp1), time + Dt);
        calculate::noise_gradient( g2, gp2, x_k + Dt*(k1+g1), phi_k + Dt*(kp1+gp1), n_x, n_phi);

        stor::phi_dw = phi_k + 0.5*(kp1 + kp2 + gp1 + gp2)*Dt;
        stor::x_dw = x_k + 0.5*(k1 + k2 + g1 + g2)*Dt;
        time += integrate::Dt;

        calculate::gradient( stor::vx, stor::phi_dt, stor::x_dw, stor::phi_dw, time);
        return 1;
    }



	double runge_kutta(double &time){

		// Store the initial positions
        phi_k = stor::phi_dw;
        x_k   = stor::x_dw;

        // Pre-calculate random numbers
        n_x = calculate::Normal();
        n_phi = calculate::Normal();

        calculate::gradient( k1, kp1, x_k, phi_k, time);
        calculate::noise_gradient( g1, gp1, x_k, phi_k, n_x, n_phi);

        // Predict a step of Dt/2 forward and calculate gradient
        calculate::gradient( k2, kp2, x_k + 0.5*Dt*(k1+g1*sqrt(2)), phi_k + 0.5*Dt*(kp1+gp1*sqrt(2)), time + 0.5*Dt);
        calculate::noise_gradient( g2, gp2, x_k + 0.5*Dt*(k1+g1*sqrt(2)), phi_k + 0.5*Dt*(kp1+gp1*sqrt(2)), n_x, n_phi);

        // Predict a step of Dt/2 forward and calculate gradient
        calculate::gradient( k3, kp3, x_k + 0.5*Dt*(k2+g2*sqrt(2)), phi_k + 0.5*Dt*(kp2+gp2*sqrt(2)), time + 0.5*Dt);
        calculate::noise_gradient( g3, gp3, x_k + 0.5*Dt*(k2+g2*sqrt(2)), phi_k + 0.5*Dt*(kp2+gp2*sqrt(2)), n_x, n_phi);

        // Predict a step of Dt forward and calculate gradient
        calculate::gradient( k4, kp4, x_k + Dt*(k3+g3*sqrt(2)), phi_k + Dt*(kp3+gp3), time + Dt);
        calculate::noise_gradient( g4, gp4, x_k + Dt*(k3+g3), phi_k + Dt*(kp3+gp3), n_x, n_phi);

        //calculate::gradient( k2, kp2, x_k + 0.5*Dt*k1, phi_k + 0.5*Dt*kp1, time + 0.5*Dt);
        //calculate::gradient( k3, kp3, x_k + 0.5*Dt*k2, phi_k + 0.5*Dt*kp2, time + 0.5*Dt);
        //calculate::gradient( k4, kp4, x_k + Dt*k3, phi_k + Dt*kp3, time + Dt);

        stor::phi_dw = phi_k + (kp1 + gp1 + 2*(kp2 +gp2 + kp3 + gp3) + kp4 + gp4)*Dt/6.0;
		stor::x_dw = x_k + (k1 + g1 + 2*(k2 + g2 + k3 + g3) + k4 + g4)*Dt/6.0;
		time += integrate::Dt;

        calculate::gradient( stor::vx, stor::phi_dt, stor::x_dw, stor::phi_dw, time);



	return 1;
	}// end of function runge_kutta

    namespace multi_dw {
        std::vector<double> x_p;
        std::vector<double> phi_p;
        std::vector<double> Kx1;
        std::vector<double> Kx2;
        std::vector<double> Kx3;
        std::vector<double> Kx4;
        std::vector<double> Kp1;
        std::vector<double> Kp2;
        std::vector<double> Kp3;
        std::vector<double> Kp4;


        void setup( int Nwires)
        {
            x_p.assign(Nwires, 0.0);
            phi_p.assign(Nwires, 0.0);
            Kx1.assign(Nwires, 0.0);
            Kp1.assign(Nwires, 0.0);
            Kx2.assign(Nwires, 0.0);
            Kp2.assign(Nwires, 0.0);
            Kx3.assign(Nwires, 0.0);
            Kp3.assign(Nwires, 0.0);
            Kx4.assign(Nwires, 0.0);
            Kp4.assign(Nwires, 0.0);
        }

        double runge_kutta(std::vector<double> &x_k, std::vector<double> &phi_k, double &time, const double dt)
        {

            calculate::gradient( Kx1, Kp1, x_k, phi_k, time);

            for( int i = 0; i < x_p.size(); i++) {
                x_p[i] = x_k[i] + 0.5*dt*Kx1[i];
                phi_p[i] = phi_k[i] + 0.5*dt*Kp1[i];
            }

            calculate::gradient( Kx2, Kp2, x_p, phi_p, time + 0.5*dt);

            for( int i = 0; i < x_p.size(); i++) {
                x_p[i] = x_k[i] + 0.5*dt*Kx2[i];
                phi_p[i] = phi_k[i] + 0.5*dt*Kp2[i];
            }

            calculate::gradient( Kx3, Kp3, x_p, phi_p, time + 0.5*dt);
            for( int i = 0; i < x_p.size(); i++) {
                x_p[i] = x_k[i] + dt*Kx3[i];
                phi_p[i] = phi_k[i] + dt*Kp3[i];
            }
            calculate::gradient( Kx4, Kp4, x_p, phi_p, time + dt);

            for( int i = 0; i < x_p.size(); i++) {
                phi_k[i] = phi_k[i] + (Kp1[i] + 2*(Kp2[i] + Kp3[i]) + Kp4[i])*dt/6.0;
                x_k[i] = x_k[i] + (Kx1[i] + 2*(Kx2[i] + Kx3[i]) + Kx4[i])*dt/6.0;
            }
            time += dt;

            //calculate::gradient( stor::vx, stor::phi_dt, stor::x_dw, stor::phi_dw, time);
            return 1;
        }// end of function runge_kutta
    }

    namespace RK45 {
        double  a[6][6] =
        { {0.0,            0.0,          0.0,           0.0,            0.0,    0.0},
          {0.25,           0.0,          0.0,           0.0,            0.0,    0.0},
          {3.0/32.0,       9.0/32.0,     0.0,           0.0,            0.0,    0.0},
          {1932.0/2197.0, -7200./2197.0, 7296.0/2197.0, 0.0,            0.0,    0.0},
          {439.0/216.0,   -8.0,          3680.0/513.0, -845.0/4104.0,   0.0,    0.0},
          {-8.0/27.0,      2.0,         -3544.0/2565.0, 1859.0/4104.0, -11./40, 0.0} };

        double b[6] = {16.0/135.0, 0.0, 6656.0/12825.0, 28561/56430.0, -9.0/50.0, 2.0/55.0};
        double bs[6] = {25.0/216.0, 0.0, 1408.0/2565.0, 2197.0/4104.0, -0.2, 0.0};
        double db[6] = { b[0]-bs[0], b[1]-bs[1], b[2]-bs[2], b[3]-bs[3], b[4]-bs[4], b[5]-bs[5]};

        double c[6] = { 0.0, 0.25, 3.0/8.0, 12.0/13.0, 1.0, 0.5};

        double kx[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        double kp[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

        double ex = 0.0, ep = 0.0;

        double runge_kutta_45( double &time, const double dt)
        {
            phi_k = stor::phi_dw;
            x_k   = stor::x_dw;

            calculate::gradient( kx[0], kp[0], x_k, phi_k, time);

            calculate::gradient( kx[1], kp[1],
                    x_k   + dt*a[1][0]*kx[0],
                    phi_k + dt*a[1][0]*kp[0],
                    time  + dt*c[1]);

            calculate::gradient( kx[2], kp[2],
                    x_k   + dt*(a[2][0]*kx[0] + a[2][1]*kx[1]),
                    phi_k + dt*(a[2][0]*kp[0] + a[2][1]*kp[1]),
                    time  + dt*c[2]);

            calculate::gradient( kx[3], kp[3],
                    x_k   + dt*(a[3][0]*kx[0] + a[3][1]*kx[1] + a[3][2]*kx[2]),
                    phi_k + dt*(a[3][0]*kp[0] + a[3][1]*kp[1] + a[3][2]*kp[2]),
                    time  + dt*c[3]);

            calculate::gradient( kx[4], kp[4],
                    x_k   + dt*(a[4][0]*kx[0] + a[4][1]*kx[1] + a[4][2]*kx[2] + a[4][3]*kx[3]),
                    phi_k + dt*(a[4][0]*kp[0] + a[4][1]*kp[1] + a[4][2]*kp[2] + a[4][3]*kp[3]),
                    time  + dt*c[4]);

            calculate::gradient( kx[5], kp[5],
                    x_k   + dt*(a[5][0]*kx[0] + a[5][1]*kx[1] + a[5][2]*kx[2] + a[5][3]*kx[3] + a[5][4]*kx[4]),
                    phi_k + dt*(a[5][0]*kp[0] + a[5][1]*kp[1] + a[5][2]*kp[2] + a[5][3]*kp[3] + a[5][4]*kp[4]),
                    time  + dt*c[5]);

            stor::x_dw   = x_k   + dt*(b[0]*kx[0] + b[1]*kx[1] + b[2]*kx[2] + b[3]*kx[3] + b[4]*kx[4] + b[5]*kx[5]);
            stor::phi_dw = phi_k + dt*(b[0]*kp[0] + b[1]*kp[1] + b[2]*kp[2] + b[3]*kp[3] + b[4]*kp[4] + b[5]*kp[5]);

            ex = dt*(db[0]*kx[0] + db[1]*kx[1] + db[2]*kx[2] + db[3]*kx[3] + db[4]*kx[4] + db[5]*kx[5]);
            ep = dt*(db[0]*kp[0] + db[1]*kp[1] + db[2]*kp[2] + db[3]*kp[3] + db[4]*kp[4] + db[5]*kp[5]);
            time += dt;

            calculate::gradient( stor::vx, stor::phi_dt, stor::x_dw, stor::phi_dw, time);
            return (ep/stor::phi_dw);
        }


    }



}// end of namespace
